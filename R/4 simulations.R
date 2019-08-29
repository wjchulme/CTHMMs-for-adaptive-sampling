
data_canon <- read_rds(path = here::here("data","data_canon.rds"))
data_wide <- read_rds(path = here::here("data","data_wide.rds"))

library('msm')




hmm_obj <- read_rds(here::here("data", "hmm_ord_2_2_3.rds"))

# choose unhealthy/alert/bad states for the model
kstar <- which(c(0, 1, 0, 1, 0, 1, 0)==1)



##

# run adaptive sampling scheme on simulated data ----------------


# sample a markov chain for 100 days, using fitted transition intensity matrix, starting at midday on day 1
MC <- sim.msm(hmm_obj$fittedQ[[1]][2:3,2:3], maxtime=100, obstime=0, mintime=0.5, start=sample(hmm_obj$n.states, size=1, prob=rep(1/hmm_obj$n.states, hmm_obj$n.states)))


set.seed(401)
#MC <- sim.msm(hmm_obj$fittedQ[[1]], maxtime=100, obstime=0, mintime=0.5, start=3)
#state_chain <- MC$states[-length(MC$states)]
#state_ttime <- MC$times[-length(MC$times)]

# simulate observation times and observations from model

adaptfunc_sim <- function(maxtime, Phi, hmodel, badstate, state_chain, state_ttime, u_0, tau, nextsample = "pass"){

  symptomsample <- Vectorize(
    function(state, hmodel){
      c(hmodel[[state]]$Anxiety$r(1),
        hmodel[[state]]$Delusions$r(1),
        hmodel[[state]]$Depression$r(1),
        hmodel[[state]]$Hallucinations$r(1)
        )
    },
    vectorize.args="state"
  )

    K <- nrow(Phi)
    Phi0 <- Phi; Phi0[badstate,] <- 0

    j <- 1
    J <- 1
    time_J <- u_0
    time <- time_J
    Y_J <- NULL
    Y <- NULL

    u_J_raw <- NULL
    u_raw <- NULL

    u_J <- NULL
    u <- NULL

    store<- NULL

    while (max(time_J) <= maxtime){

      truestate <- state_chain[findInterval(time_J, state_ttime)]

      Y_J <- matrix(symptomsample(truestate, hmodel), nrow=1)
      Y <- rbind(Y, Y_J)


      p_1J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=J)

      p_J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=1)


      ## based on snapshot risk
      u_J <- 3/24 #must sample at least three hours after last sample

      if(is.numeric(nextsample)){
        #sample deterministically
        u_J_raw <- nextsample
      } else {
        u_J_raw <- select_next_time(p_1J, Phi, tau, badstate)
      }

      u_J <- adjust_next_time(time_J, u_J_raw) - time_J


      u_raw <- c(u_raw, u_J_raw)
      u <- c(u, u_J)

      store_J <-
        tibble(
          tau = tau,
          J = J,
          j = list(j),
          time_J = time_J,
          time = list(time),
          y_J = list(Y_J),
          y = list(Y),
          u_J_raw = u_J_raw,
          u_J = u_J,

          p_J = list(p_J),

          p_smooth = pmap(list(time, y, J), function(time, y, J){stateprobMV_FB(Phi, hmodel.in, time, y, J)}),
          #p_filter = accumulate(p_smooth, ~bind_rows(.x, slice(.y, n()))),
          p_1J = map(p_smooth, function(p_smooth){unlist(slice(p_smooth[,-1], n()))}),

          #u = map(time_J, ~seq(0.2, (maxtime+0.2) - ., +0.2)),

          # PS = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
          #     mat<-futurestateprob(p_1J, Phi, u)
          #     rownames(mat)<-names(p_1J)
          #     mutate(as_tibble(t(mat)), time=time_J+u)
          #   }),
          # PP = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
          #     mat<-futurepassprob(p_1J, Phi, u)
          #     rownames(mat)<-names(p_1J)
          #     mutate(as_tibble(t(mat)), time=time_J+u)
          #   }),


        ) %>% select(-p_smooth, -time, -j, -Y)

      store <- bind_rows(store, store_J)

      J <- J+1
      j <- seq_len(J)
      time_J <- time_J + u_J
      time <- c(time, time_J)

    }


    store

}

hmodel.in <-
  hmm_obj$problist[[1]] %>%
  map_depth(1, ~set_names(., c("Anxiety", "Depression", "Delusions", "Hallucinations"))) %>%
  map_depth(2, ~msm::hmmCat(prob = .)) %>% map_depth(1, ~lift_dl(msm::hmmMV)(.))



## run this for:
#statestate0 = 1; badstate0=2
#statestate0 = 3; badstate0=4
statestate0 = 5; badstate0=6



design <-
  cross_df(list(
    Phi = list(list(hmm_obj$fittedQ[[1]])),
    hmodel.in = list(hmodel.in),
    startstate = startstate0,
    badstate = list(c(badstate0)),
    repID = 1:500,
    maxtime = 100,
  )) %>%
  mutate(Phi=map(Phi, ~.[[1]]))

design2 <- design %>%
  mutate(
    MC = pmap(list(Phi, maxtime, startstate), function(Phi, maxtime, startstate){sim.msm(Phi, maxtime=maxtime, obstime=0, mintime=9/24, start=startstate)}),
    state_chain = map(MC, ~ .$states[-length(.$states)]),
    state_ttime = map(MC, ~ .$times[-length(.$times)]),
    badstate_ttime = pmap(list(state_chain, state_ttime, badstate), ~..2[..1 %in% ..3]),
  )

design3 <- purrr::map_dfr(c(0, 1/32, 1/16, 1/8, 1/4, 1/2) %>% set_names(.,.) , ~design2, .id="tau") %>%
  mutate(tau=as.numeric(tau))



## send this to iCSF as will tkae a long time to run
library('furrr')
parallel::detectCores()
n_cores <- 10

clusters <- parallel::makeCluster(n_cores)
plan(cluster, workers=clusters)

tictoc::tic()
simresults <- design3 %>%
  mutate(
     adp = future_pmap(list(Phi, hmodel.in, badstate, tau, maxtime, state_chain, state_ttime),
               function(Phi, hmodel.in, badstate, tau, maxtime, state_chain, state_ttime){
                 adaptfunc_sim(maxtime=maxtime, Phi=Phi, hmodel=hmodel.in, badstate=badstate, state_chain=state_chain, state_ttime=state_ttime, u_0=12/24, tau=tau, nextsample = "pass")
                 }
     ),
 )
tictoc::toc()

parallel::stopCluster(clusters)
plan(sequential)

write_rds(simresults, path = here::here("data", "simulations", str_c("simresults_223_ss",startstate0,"_bs",badstate0,".rds")), compress="xz")


# run adaptive sampling scheme on real data ----------------

adaptfunc_real <- function(Phi, hmodel, Y_complete, time_complete, badstate, tau, nextsample = "pass"){

    maxtime = max(time_complete)
    K <- nrow(Phi)
    Phi0 <- Phi; Phi0[badstate,] <- 0

    J <- 4
    time_J <- time_complete[J]
    time <- time_complete[seq_len(J)]
    Y_J <- NULL
    Y <- NULL

    u_J_raw <- NULL
    u_raw <- NULL

    u_J <- NULL
    u <- NULL

    store<- NULL

    while (max(time_J) <= maxtime){

      #Y_J <- Y_complete[J,]
      Y <- Y_complete[seq_len(J),]

      p_1J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=J)

      #p_J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=2)

      #u_J_raw <- select_next_time_incr(p_1J, Phi, tau, badstate, increment=30/24/60, max=3.1)
      u_J_raw <- select_next_time(p_1J, Phi, tau, badstate, tol=0.0001)
      u_raw <- c(u_raw, u_J_raw)

      if(time_J+u_J_raw >max(time_complete)){break}

      time_Jplus1 <- findclosest(time_J+u_J_raw, time_complete)


      u_J <- time_Jplus1 - time_J
      u <- c(u, u_J)

      store_J <-
        tibble(
          tau = tau,
          J = J,
          time_J = time_J,
          time = list(time),
          y_J = list(Y_J),
          y = list(Y),
          u_J_raw = u_J_raw,
          u_J = u_J,

          #p_J = list(p_J),

          #p_smooth = pmap(list(time, y, J), function(time, y, J){stateprobMV_FB(Phi, hmodel.in, time, y, J)}),
          #p_filter = accumulate(p_smooth, ~bind_rows(.x, slice(.y, n()))),
          p_1J = list(p_1J), #map(p_smooth, function(p_smooth){unlist(slice(p_smooth[,-1], n()))}),

          #u = map(time_J, ~seq(0.2, (maxtime+0.2) - ., +0.2)),

          # PS = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
          #     mat<-futurestateprob(p_1J, Phi, u)
          #     rownames(mat)<-names(p_1J)
          #     mutate(as_tibble(t(mat)), time=time_J+u)
          #   }),
          # PP = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
          #     mat<-futurepassprob(p_1J, Phi, u)
          #     rownames(mat)<-names(p_1J)
          #     mutate(as_tibble(t(mat)), time=time_J+u)
          #   }),


        ) #%>% select(-p_smooth)

      store <- bind_rows(store, store_J)

      J <- which(abs(time_complete-time_Jplus1)<0.001)
      time_J <- time_complete[J]
      time <- time_complete[seq_len(J)]

    }

    store

}

hmodel.in <-
  hmm_obj$problist[[1]] %>%
  map_depth(1, ~set_names(., c("Anxiety", "Depression", "Delusions", "Hallucinations"))) %>%
  map_depth(2, ~msm::hmmCat(prob = .)) %>% map_depth(1, ~lift_dl(msm::hmmMV)(.))


design <-
  cross_df(list(
    Phi = list(list(hmm_obj$fittedQ[[1]])),
    hmodel.in = list(hmodel.in),
    #tau = c(0, 0.1, 0.3, 0.5),
    badstate = list(list(kstar)),
    ptid = list(unique(data_wide$ptid))
  )) %>%
  mutate(
    Phi=map(Phi, ~.[[1]]),
    badstate=map(badstate, ~.[[1]])
  )


design2 <- purrr::map_dfr(c(0, 1/32, 1/16, 1/8, 1/4, 1/2) %>% set_names(.,.) , ~design, .id="tau") %>%
  mutate(tau=as.numeric(tau))


library('furrr')
parallel::detectCores()
n_cores <- 10

clusters <- parallel::makeCluster(n_cores)
plan(cluster, workers=clusters)

tictoc::tic()
realresults <- design2 %>%
  mutate(

    adp = pmap(list(Phi, hmodel.in, badstate, tau, ptid),

                 possibly(
                    function(Phi, hmodel.in, badstate, tau, ptid2){
                      dat <- filter(data_wide, ptid==ptid2)
                      adaptfunc_real(Phi=Phi, hmodel=hmodel.in, Y_complete=dat$responseMV, time_complete = dat$fup_day, badstate=badstate, tau=tau, nextsample = "pass")
                    },
                    list(NA)
                 )
    ),
  )
tictoc::toc()

parallel::stopCluster(clusters)
plan(sequential)

write_rds(realresults, path=here::here("data", "simulations", "realresults_223.rds"), compress="gz")



