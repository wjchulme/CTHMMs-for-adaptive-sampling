library('tidyverse')
library('msm')




# number of states
K <- 3

# rate parameters
Theta_Anx <-
  rbind(
    c(0.665994977, 0.142934090, 0.134067351, 0.032646402, 0.002760432, 0.003723969, 0.017872780),
    c(9.997197e-01, 2.314477e-04, 4.704221e-05, 1.800715e-06, 1.226944e-08, 0.000000e+00, 0.000000e+00),
    c(0.407222161, 0.176092029, 0.237077077, 0.084087767, 0.008329331, 0.011722988, 0.075468648)
  )

Theta_Del <-
  rbind(
    c(0.02934004, 0.03009703, 0.23826052, 0.18248337, 0.15933683, 0.30072651, 0.05975569),
    c(0.11885533, 0.07890069, 0.37320357, 0.17446239, 0.11159712, 0.13127620, 0.01170470),
    c(0.12025091, 0.07945592, 0.37399733, 0.17395682, 0.11093152, 0.12991492, 0.01149258)
  )

Theta_Dep <-
  rbind(
    c(1.217296e-08, 5.077457e-06, 8.032691e-05, 1.563228e-04, 6.463392e-04, 9.255951e-02, 9.065524e-01),
    c(9.958148e-01, 3.560538e-03, 5.707208e-04, 3.679823e-05, 1.407788e-05, 3.055733e-06, 0.000000e+00),
    c(0.6297103812, 0.1919219126, 0.1196058648, 0.0254077475, 0.0193742179, 0.0139489226, 0.0000309534)
  )


Theta_Hal <-
  rbind(
    c(0.26218142, 0.13607401, 0.11547661, 0.41639301, 0.02306948, 0.02759212, 0.01921335),
    c(0.34464100, 0.14695777, 0.11529528, 0.34979831, 0.01546282, 0.01731972, 0.01052510),
    c(0.04698227, 0.05049253, 0.06026716, 0.51175076, 0.06891337, 0.11057198, 0.15102194)
  )



Theta <- list(Theta_Anx, Theta_Del, Theta_Dep, Theta_Hal)



# transition intensities
Phi <-
  rbind(
    c(-0.52255092,  0.24877507,  0.27377585),
    c( 0.08138638, -0.1805877,  0.09920132),
    c( 0.07028803,  0.1020642, -0.17235223)
  )

# state starting probabilities
pi_1 <- c(1,0,0)








# compose emission distribution parameters in form needed for msm

Theta_arr <- array(c(Theta_Anx, Theta_Del, Theta_Dep, Theta_Hal),
                   dim=c(K, 7L, 4L),
                   dimnames=
                     list(
                       state=paste0("S", 1:K),
                       value=1:7,
                       item=c("Anxiety", "Delusions", "Depression", "Hallucinations")
                     )
                   )

Theta_list <- array_tree(Theta_arr, margin=c(1,3))

hmodel.in <-
  list(
    lift_dl(hmmMV)(
      map(Theta_list[[1]], ~hmmCat(prob=., basecat=1))
    ),
    lift_dl(hmmMV)(
      map(Theta_list[[2]], ~hmmCat(prob=., basecat=1))
    ),
    lift_dl(hmmMV)(
      map(Theta_list[[3]], ~hmmCat(prob=., basecat=1))
    )
  )





# simulate states and responses for above model

obs.scheme <- tibble(time = 1:400)

set.seed(420)
simobj <-
  simmulti.msm(
    data = obs.scheme,
    start = 1,
    qmatrix = Phi#,
    #hmodel = hmodel.in
  ) %>%
  select(-subject, -keep) %>%
  mutate(
    Anxiety = map_dbl(state, ~sample(1:7, size=1, replace=TRUE, prob = Theta_Anx[.,])),
    Delusions = map_dbl(state, ~sample(1:7, size=1, replace=TRUE, prob = Theta_Del[.,])),
    Depression = map_dbl(state, ~sample(1:7, size=1, replace=TRUE, prob = Theta_Dep[.,])),
    Hallucinations = map_dbl(state, ~sample(1:7, size=1, replace=TRUE, prob = Theta_Hal[.,]))
  )

simobj$obsMV <- cbind(simobj$Anxiety, simobj$Delusions, simobj$Depression, simobj$Hallucinations)

# create hmm model object
hmmobj <-
  msm(
    data = simobj,
    formula = obsMV ~ time,
    obstype = 1,
    #initprobs = pi_1,
    qmatrix = Phi,
    hmodel = hmodel.in,
    fixedpars = TRUE
  )


stateprobMV <- function(Phi, hmodel, time, y, w){

  if(!is.matrix(y)) {
    y <- matrix(y, nrow=1)
  }
  if((length(time) != nrow(y))) stop("dimensions of time and y not compatible")
  if(w<1 | w>length(time)) stop("w must be between 1 and the number of observations available")

  J <- length(time)
  time <- time[seq_len(w)+(J-w)]
  y <- y[seq_len(w)+(J-w), ]


  if(w==1) {

    #This is a mixture model since we only have one time-point (one observation)

    S <- length(hmodel)
    weights <- rep(NA, S)
    for(s in seq_len(S)){
      weights[s] <- (hmodel.in[[s]]$Anxiety$pars[3:9][y[1]] +
                     hmodel.in[[s]]$Delusions$pars[3:9][y[2]] +
                     hmodel.in[[s]]$Depression$pars[3:9][y[3]] +
                     hmodel.in[[s]]$Hallucinations$pars[3:9][y[4]])/4
    }
    p_tJ <- weights/sum(weights)
  } else {

      #we create a msm object for convenience so we can call the viterbi.msm
      #which, confusingly, also runs the forward-backward algorithm, which
      #calculates the probability of each state at each time point conditional on all the data
      #we take the probabilities for the last time-point (equivalent to the forward algorithm!)
      obj <- msm(
        formula = y ~ time,
        #initprobs = initprobs,
        est.initprobs = TRUE,
        qmatrix = Phi,
        hmodel = hmodel,
        fixedpars = TRUE
      )

      p_t <- viterbi.msm(obj)$pstate
      p_tJ <- p_t[nrow(p_t),]
  }
  names(p_tJ) <- paste0("S",seq_len(nrow(Phi)))
  p_tJ
}



stateprobMV_FB <- function(Phi, hmodel.in, time, y, w){

  if(!is.matrix(y)) {
    y <- matrix(y, nrow=1)
  }
  if((length(time) != nrow(y))) stop("dimensions of time and y not compatible")
  if(w<1 | w>length(time)) stop("w must be between 1 and the number of observations available")

  J <- length(time)
  time <- time[seq_len(w)+(J-w)]
  y <- y[seq_len(w)+(J-w), ]


  if(w==1) {

    #This is a mixture model since we only have one time-point (one observation)

    S <- length(hmodel.in)
    weights <- rep(NA, S)
    for(s in seq_len(S)){
      weights[s] <- (hmodel.in[[s]]$Anxiety$pars[3:9][y[1]] +
                     hmodel.in[[s]]$Delusions$pars[3:9][y[2]] +
                     hmodel.in[[s]]$Depression$pars[3:9][y[3]] +
                     hmodel.in[[s]]$Hallucinations$pars[3:9][y[4]])/4
    }
    p_t <- t(weights/sum(weights))
    colnames(p_t) <- paste0("S",seq_len(nrow(Phi)))
    bind_cols(tibble(time), as_tibble(p_t))
  } else {

      #we create a msm object for convenience so we can call the viterbi.msm
      #which, confusingly, also runs the forward-backward algorithm, which
      #calculates the probability of each state at each time point conditional on all the data
      #we take the probabilities for the last time-point (equivalent to the forward algorithm!)
      obj <- msm(
        formula = y ~ time,
        #initprobs = initprobs,
        est.initprobs = TRUE,
        qmatrix = Phi,
        hmodel = hmodel.in,
        fixedpars = TRUE
      )


    vit <- viterbi.msm(obj)
    time <- vit$time
    p_t <- vit$pstate
    colnames(p_t) <- paste0("S",seq_len(nrow(Phi)))
    bind_cols(tibble(time), as_tibble(p_t))
  }

}


futurestateprob <- Vectorize(function(p, Phi, u){
  p %*% expm::expm(u*Phi)
}, vectorize.args="u")

futurepassprob<-Vectorize(function(p, Phi, u){
  ppass.msm(qmatrix = Phi, start=p, tot=u)
}, vectorize.args="u")

firstpassagetime <- function(p, Phi, tostate){
  efpt.msm(qmatrix = Phi, start=p, tostate=tostate)
}



############# adaptive sampl

## simulate state sequence from model

set.seed(401)
obs.scheme <- tibble(time = 1)
MC <- sim.msm(Phi, maxtime=365, obstime=0, mintime=0, start=sample(K, size=1, prob=c(0, 0, 1)))

state_chain <- MC$states
state_ttime <- MC$times

# simulate observation times and observations from model

adaptfunc <- function(maxtime, Phi, hmodel, state_chain, state_ttime, ustar_0, state_star, nextsample=c("pass","prob")){

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
  Phi0 <- Phi; Phi0[state_star,] <- 0

  j <- 1
  J <- 1
  time_J <- ustar_0
  time <- time_J
  Y_J <- NULL
  Y <- NULL
  ustar_J <- NULL
  ustar <- NULL

  store<- NULL

  while (max(time_J) <= maxtime){

    truestate <- state_chain[findInterval(time_J, state_ttime)]

    Y_J <- matrix(symptomsample(truestate, hmodel), nrow=1)
    Y <- rbind(Y, Y_J)


    p_1J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=J)

    p_J <- stateprobMV(Phi=Phi, hmodel=hmodel, time = time, y=Y, w=1)


    ## based on snapshot risk
    ustar_J <- 2/24 #must sample at least two hours after last sample

    if(is.numeric(nextsample)){
      #sample deterministically
      ustar_J<-nextsample
    }


    if(nextsample=="prob"){

      p_1J_ustar <- futurestateprob(p_1J, Phi, ustar_J)
      while(
        ( sum(p_1J_ustar[state_star]) < 0.1 |           #must not sample when combined risk of states 3, 6, or 9 does not exceed p=0.2
          !between((time_J+ustar_J)%%1, 9/24, 21/24)) &  #must not sample outside of 9am to 9pm
        ustar_J < 3                                  #must sample before 72 hours have passed
      ){
        ustar_J <- ustar_J+(10/24/60) # consider in 10 minute increments
        p_1J_ustar <- futurestateprob(p_1J, Phi, ustar_J)
      }
    }

    if(nextsample=="pass"){

      #based on cumulative risk - TODO - need to figure out how to sum these non-independent probabilities
      p_1J_ustar <- futurepassprob(p_1J, Phi, ustar_J)

      while(
        ( sum(p_1J_ustar[state_star]) < 0.1 |           #must sample when combined risk of states 3, 6, or 9 exceeds p=0.2
          !between((time_J+ustar_J)%%1, 9/24, 21/24)) &  #must not sample outside of 9am to 9pm
        ustar_J < 3                                  #must sample before 72 hours have passed
      ){
        ustar_J <- ustar_J+(10/24/60) # consider in 10 minute increments
        p_1J_ustar <- futurepassprob(p_1J, Phi0, ustar_J)
      }
    }

    ustar <- c(ustar, ustar_J)

    store_J <-
      tibble(
        J = J,
        j = list(j),
        time_J = time_J,
        time = list(time),
        y_J = list(Y_J),
        y = list(Y),
        ustar_J = ustar_J,
        ustar = list(ustar),
        state = list(seq_len(K)),
        p_J = list(p_J),

        p_smooth = pmap(list(time, y, J), function(time, y, J){stateprobMV_FB(Phi, hmodel.in, time, y, J)}),
        p_filter = accumulate(p_smooth, ~bind_rows(.x, slice(.y, n()))),
        p_1J = map(p_smooth, function(p_smooth){unlist(slice(p_smooth[,-1], n()))}),

        u = map(time_J, ~seq(0.2, (maxtime+0.2) - ., +0.2)),

        PS = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
            mat<-futurestateprob(p_1J, Phi, u)
            rownames(mat)<-names(p_1J)
            mutate(as_tibble(t(mat)), time=time_J+u)
          }),
        PP = pmap(list(p_1J, u, time_J), function(p_1J, u, time_J){
            mat<-futurepassprob(p_1J, Phi, u)
            rownames(mat)<-names(p_1J)
            mutate(as_tibble(t(mat)), time=time_J+u)
          }),


      )

    store <- bind_rows(store, store_J)

    J <- J+1
    j <- seq_len(J)
    time_J <- time_J + ustar_J
    time <- c(time, time_J)

  }


  store

}


state_probs_adp <- adaptfunc(100, Phi, hmodel.in, state_chain, state_ttime, ustar_0=10/24, state_star=c(1), nextsample="pass")






#########################################







state_probs<-state_probs_adp

  png(file=".\\figures\\gifs\\temp%03d.png", width=600, height=600)
  indexJ <- 0
  nexttime_J <- 0
    while (nexttime_J <=100){

    indexJ <- indexJ+1
    currentime_J <- state_probs$time_J[[indexJ]]
    nexttime_J <- state_probs$time_J[[indexJ]]+state_probs$ustar_J[[indexJ]]

    plot_obs <- state_probs %>%
    filter(J==indexJ) %>%
    mutate(obs=map2(y, time, ~bind_cols(tibble(time=.y), as_tibble(.x) %>% set_names(c("Anxiety","Delusions","Depression","Hallucinations"))))) %>%
    select(obs) %>%
    unnest() %>%
    gather("item", "obs", -time) %>%
    ggplot() +
    facet_grid(rows=vars(item))+
    geom_point(aes(x=time, y=obs, colour=factor(obs, 0:7)))+
    labs(colour="", x='time', y='response')+
    scale_colour_manual(values=set_names(c("lightgrey",c(viridisLite::viridis(n=7))), 0:7), na.value="transparent", breaks=0:7, labels=c("no response",1:7))+
    scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
    scale_y_continuous(breaks=c(1,4,7), limits=c(0,8))+
    theme_bw()+
    theme(
      legend.position="bottom",legend.direction="horizontal",
      #strip.text = element_blank(),
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    guides(colour=guide_legend(nrow=1))

#
#     plot_filter <- state_probs %>%
#     filter(J==indexJ) %>%
#     mutate(
#       p_filter_proj = map2(p_filter, PS, bind_rows)
#     ) %>%
#     select(p_filter_proj) %>%
#     unnest() %>%
#     gather("state", "prob", -time) %>%
#     ggplot() +
#     geom_area(aes(x=time, y=prob, fill=as.factor(state)), colour='transparent', position='stack')+
#     geom_vline(aes(xintercept=currentime_J), linetype='dashed')+
#     labs(y='PS', fill='state', y=NULL)+
#     scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
#     #scale_fill_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
#     scale_fill_manual(values=c("#e6550d", "#fdae6b", "#fee6ce"))+
#     scale_linetype(guide=FALSE)+
#     theme_bw()+
#     theme(legend.position="none")
#


    plot_smooth <- state_probs %>%
    filter(J==indexJ) %>%
    select(p_smooth) %>%
    unnest() %>%
    gather("state", "prob", -time) %>%
    ggplot() +
    geom_area(aes(x=time, y=prob, fill=as.factor(state)), colour='transparent', position='stack')+
    geom_vline(aes(xintercept=currentime_J), linetype='dashed')+
    labs(y='Pr(z_t|y_1:T)', fill='state', y=NULL)+
    scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
    #scale_fill_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
    scale_fill_manual(values=c("#e6550d", "#fdae6b", "#fee6ce"))+
    #scale_y_continuous(limits=c(0,0.40))+
    scale_linetype(guide=FALSE)+
    theme_bw()



    plot_pass <- state_probs %>%
    filter(J==indexJ) %>%
    select(PP) %>%
    unnest() %>%
    gather("state", "prob", -time) %>%
    ggplot() +
    geom_line(aes(x=time, y=prob, colour=as.factor(state)))+
    geom_vline(aes(xintercept=currentime_J), linetype='dashed')+
    labs(y='PP', colour='state', y=NULL)+
    scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
    #scale_colour_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
    scale_colour_manual(values=c("#e6550d", "#fdae6b", "#fee6ce"))+
    #scale_y_continuous(limits=c(0,0.40))+
    scale_linetype(guide=FALSE)+
    theme_bw()


    print(cowplot::plot_grid(plot_obs, plot_smooth, plot_pass, ncol=1, align="v", axis="lr", rel_heights=c(2,1,1)))


  }

dev.off()
# convert pngs to one gif using ImageMagick
#system("magick convert -delay 20 *.png example_2_reduced.gif", intern=TRUE)

system('"C:\\Program Files\\ImageMagick-7.0.8-Q16\\magick.exe" -delay 20 figures/gifs/*.png "figures/gifs/adp_sampling_v8.gif"')
# cleaning up
file.remove(list.files(path=here::here("figures","gifs"), pattern=".png", full.names=TRUE))








library('gganimate')

  plot_obs <- state_probs %>%
    mutate(obs=map2(y, time, ~bind_cols(tibble(time=.y), as_tibble(.x) %>% set_names(c("Anxiety","Delusions","Depression","Hallucinations"))))) %>%
    select(obs, time_J, J) %>%
    unnest() %>%
    gather("item", "obs", -time, -time_J, -J) %>%
    ggplot() +
    facet_grid(rows=vars(item))+
    geom_point(aes(x=time, y=obs, colour=factor(obs, 0:7)))+
    transition_time(time_J)+
    view_static()+ease_aes('linear')+
    labs(colour="", x='time', y='response')+
    scale_colour_manual(values=set_names(c("lightgrey",c(viridisLite::viridis(n=7))), 0:7), na.value="transparent", breaks=0:7, labels=c("no response",1:7))+
    scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
    scale_y_continuous(breaks=c(1,4,7), limits=c(0,8))+
    theme_bw()+
    theme(
      legend.position="bottom",legend.direction="horizontal",
      #strip.text = element_blank(),
      strip.text.y = element_text(angle = 0),
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    guides(colour=guide_legend(nrow=1))



anim_data <-
  state_probs %>%
  select(p_smooth, time_J, J) %>%
  unnest(p_smooth) %>%
  gather("state", "prob", -time, -time_J, -J) %>%
  arrange(J, time, state)


plotanim_smooth <- anim_data %>%
  ggplot() +
  geom_area(aes(x=time, y=prob, fill=as.factor(state)), colour='transparent', position='stack')+
  geom_vline(aes(xintercept=time_J), linetype='dashed')+
  transition_time(time_J)+
  view_static()+ease_aes('linear')+
  labs(y='Pr(z_t|Y_1:t)', fill='state', y=NULL)+
  scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+
  scale_fill_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
  scale_linetype(guide=FALSE)+
  theme_bw()+
  theme(legend.position="none")


plotanim_filter <- anim_data %>%
  ggplot() +
  geom_area(aes(x=time, y=prob, fill=as.factor(state)), colour='transparent', position='stack')+
  geom_vline(aes(xintercept=time_J), linetype='dashed')+
  transition_time(time_J) +
  view_static() + ease_aes('quadratic-in')+
  labs(y='Pr(z_t|Y_1:t)', fill='state', y=NULL)+
  scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100))+
  scale_fill_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
  scale_linetype(guide=FALSE)+
  theme_bw()+
  theme(legend.position="none")



  plot_pass <- state_probs %>%
    select(J, time_J, PP) %>%
    unnest() %>%
    gather("state", "prob", -time, -time_J, -J) %>%
    ggplot() +
    geom_line(aes(x=time, y=prob, colour=as.factor(state)))+
    geom_vline(aes(xintercept=time_J), linetype='dashed')+
    transition_time(time_J) +
    view_static() + #ease_aes('linear')+
    labs(y='PP', colour='state', y=NULL)+
    scale_x_continuous(breaks=seq(0,100,7), limits=c(0,100))+
    #scale_colour_manual(values=c("#e5f5e0", "#a1d99b", "#31a354", "#fee6ce", "#fdae6b", "#e6550d", "#efedf5", "#bcbddc", "#756bb1"))+
    scale_colour_manual(values=c("#e6550d", "#fdae6b", "#fee6ce"))+
    #scale_y_continuous(limits=c(0,0.40))+
    scale_linetype(guide=FALSE)+
    theme_bw()
