
source(here::here("R", "0 preliminaries.R"))

## modified version of msm package to include truncated normal + ordinal regression response distributions
remove.packages("msm")
install.packages("C:\\Users\\William\\OneDrive\\Work\\packages\\msm_1.6.7.tar.gz", repos = NULL, type="source")
library('msm')




# prepare data for use with msm ----------------------------------------------------------------------


data_canon <- read_rds(path = here::here("data","data_canon.rds"))

data_wide <-
  data_canon %>% select(ptid, fup_days, item = itemgrp, response = score) %>%
  group_by(ptid, fup_days) %>%
  spread(item, response) %>%
  ungroup() %>%
  arrange(ptid, fup_days) %>%
  ungroup()

data_wide$responseMV <- cbind(Anxiety = data_wide$Anxiety,
                              Delusions = data_wide$Delusions,
                              Hallucinations = data_wide$Hallucinations,
                              Depression = data_wide$Depression
                              )




###################################################################################################
# model configuration -------- --------------------------------------------------------------------
###################################################################################################
#
#. roughly means a function
#_ roughly means a subscript

# one row per model configuration
#
configuration0 <-
  tibble(
    n.clusters = 1:4 #number of clusters
  ) %>%
  mutate(
    n.states_cluster = map(n.clusters, ~arrangements::combinations(1:9, .x, replace=TRUE, layout="list")) #number of states in each cluster
  ) %>%
  unnest(n.states_cluster) %>%
  mutate(
    n.states = map_dbl(n.states_cluster, sum) #number of states overall
    ) %>%
  mutate(
    chr.n.states_cluster = map_chr(n.states_cluster, ~str_c(as.character(.), collapse="_")), #number of states per cluster (unique character ID for model)
    index_state = map(n.states, ~str_c("s",seq_len(.x))), #index of the state
    index_cluster_state = map(n.states_cluster, ~rep.int(seq_along(.x), times=.x)), #the cluster of the state
  ) %>%
  select(n.clusters, n.states, n.states_cluster, chr.n.states_cluster, index_state, index_cluster_state) %>%
  arrange(n.clusters, n.states) %>%
  mutate(
    name_channel = list(c("Anxiety",
                          "Delusions",
                          "Hallucinations",
                          "Depression"
                          )
                        ), # name of the channel/item
    index_channel = map(name_channel, seq_along),
    n.channels = map_dbl(name_channel, length),
    formula = map(name_channel, ~as.formula(str_c("cbind(",str_c(., collapse = ","), ")~ fup_days"))),

    n.cats = 7,
    n.cats_channel = map(name_channel, ~setNames(rep(7, length(.)), .)), # the support (number of categories) of the channel
    index_cats_channel = map(name_channel, ~setNames(rep(list(1:7), length(.)), .)) #the index of the categories for each channel
  ) %>%
  select(chr.n.states_cluster, formula, n.channels, n.states, n.clusters, n.states_cluster, n.cats, index_state, index_cluster_state, index_channel, name_channel, n.cats_channel, index_cats_channel) %>%
  filter(
    #consider between 2 and 3 states per cluster and between 1 to 3 clusters in total
    n.clusters %in% c(1,2,3),
    map_lgl(n.states_cluster, ~all(. %in% c(2,3)))
  )

write_rds(configuration0, path = here::here("data", str_c("configuration0.rds")))

#read_rds(path = here::here("data", str_c("configuration0 categorical.rds")))






###################################################################################################
# prepare categorical models  --------------------------------------------------------------------
###################################################################################################


seed <- 123
set.seed(seed)



# create random/fixed starting values and other model inputs for msm::msm() call
configuration_cat <-
  configuration0 %>%
  mutate(

    emissiondist = list(partial(msm::hmmCat, basecat=1)),
    emissiondist_channel = map(name_channel, ~setNames(rep(list(partial(msm::hmmCat, basecat=1)), length(.)), .)), #emission function of the channel

  #non-random starting values
    #init.state.probs=map(n.states, ~rep(1/.x, .x)), #starting values evenly spread
    init.trans.Q = map(n.states_cluster, ~make.init.fixed.Q(., out.rate=0.1)), #transition intensity
    init.trans.P = map(init.trans.Q, ~expm::expm(.)), #transition probability
    #init.params.fixed_channel_state = pmap(list(index_state, name_channel, n.cats),
    #                                  function(index_state, name_channel, n.cats){
    #                                    array(
    #                                      rep(1/n.cats, n.cats*length(name_channel)*length(index_state)),
    #                                      dim = c(n.cats, length(name_channel), length(index_state)),
    #                                      dimnames = list('prob'=seq_len(n.cats), 'channel'=name_channel, 'state'=index_state)
    #                                    ) %>%
    #                                   array_tree(margin = c(3,2))
    #                                    }
    #                                  ),

  #random starting values
    init.state.probs = map(n.states, ~dirichlet::rdirichlet(1, rep(1,.x)) %>% as.vector()), #random starting values
    #init.trans.P = map(n.states_cluster, ~map(.x, function(states){(dirichlet::rdirichlet(states, rep(1,states))+diag(states)) %>% {./rowSums(.)}}) %>% Matrix::bdiag() %>% as.matrix()), #random transition probability
    #init.trans.Q = map(init.trans.P, ~expm::logm(., method='Eigen')), #transition intensity
    init.params.random_channel_state = pmap(list(index_state, name_channel, n.cats),
                                     function(index_state, name_channel, n.cats){
                                       dirichlet::rdirichlet(length(name_channel)*length(index_state), rep(1, n.cats)) %>%
                                       t() %>% unlist() %>%
                                       array(
                                         dim = c(n.cats, length(name_channel), length(index_state)),
                                         dimnames = list('prob'=seq_len(n.cats), 'channel'=name_channel, 'state'=index_state)
                                       ) %>% array_tree(margin = c(3,2))
                                       }
                                     ),




    hmodel.in = pmap(list(emissiondist, init.params.random_channel_state),
                  function(emissiondist, params_channel_state){
                    params_channel_state %>% map_depth(2, ~emissiondist(prob = .)) %>% map_depth(1, ~lift_dl(msm::hmmMV)(.))
                    }
                  ),

    n.pars = pmap_dbl(list(n.cats, n.channels, n.states_cluster, n.states, n.clusters),
                 function(n.cats, n.channels, n.states_cluster, n.states, n.clusters){
                   n.qpars <- sum(n.states_cluster*(n.states_cluster-1))
                   n.catpars <- n.states*n.channels*(n.cats-1)
                   n.estinitprobs <- (n.states-1)
                   n.qpars+n.catpars+n.estinitprobs
                   }
              ),

)


# fit categorical hmms from list of inputs in configuration ---------------------------------------------

## send this to iCSF as this will take a LONG time to run

library('furrr')
parallel::detectCores()
n_cores <- 10
clusters <- parallel::makeCluster(n_cores)
plan(cluster, workers=clusters)

tictoc::tic()
hmmslist_cat0 <-
  configuration_cat %>%
  mutate(
    hmmfit = future_pmap(
      list(chr.n.states_cluster, init.state.probs, init.trans.Q, hmodel.in, n.pars),
      function(chr.n.states_cluster, init.state.probs, init.trans.Q, hmodel.in, n.pars)
      {
        start_time <- Sys.time()
        mod<-msm::msm(
          formula = responseMV ~ fup_days,
          subject = ptid,
          data = data_wide,
          initprobs = init.state.probs,
          est.initprobs = TRUE,
          qmatrix = init.trans.Q,
          gen.inits = FALSE,
          hmodel = hmodel.in,
          obstype = 1,
          opt.method = "optim",
          control = list(fnscale=10000,
                         maxit=10000, #default=500 for Nelder mead
                         reltol=1e-10, #default~=1e-08
                         ndeps = rep(1e-05, n.pars) #default=1e-3
                         ),
           method="Nelder-Mead" #default="BFGS"
        )

      }#,.progress = TRUE
    ),
  )
tictoc::toc()
parallel::stopCluster(clusters)
plan(sequential)


write_rds(hmmslist_cat0, path = here::here("data", str_c("hmmslist0 categorical grp seed is ", seed,".rds")), compress = "gz")





###################################################################################################
# prepare ordinal models  --------------------------------------------------------------------
###################################################################################################


seed <- 123
set.seed(seed)


configuration_ord <-
  configuration0 %>%
  mutate(

    trunc = list(5),
    emissiondist = map(trunc, ~(partial(msm::hmmClmTNorm7, trunc = .x))),
    emissiondist_channel = map2(name_channel, trunc, ~setNames(rep(list(partial(msm::hmmClmTNorm7, trunc = .y)), length(.x)), .x)), #emission function of the channel, if this differs from channel to channel

  #non-random starting values
    #init.state.probs=map(n.states, ~rep(1/.x, .x)), #starting values evenly spread
    init.trans.Q = map(n.states_cluster, ~make.init.fixed.Q(., out.rate=0.1)), #transition intensity
    init.trans.P = map(init.trans.Q, ~expm::expm(.)), #transition probability
    #init.params.fixed_channel_state = pmap(list(index_state, name_channel, n.cats, trunc),
    #                                  function(index_state, name_channel, n.cats, trunc){
    #                                    array(
    #                                      rep(1/n.cats, n.cats*length(name_channel)*length(index_state)),
    #                                      dim = c(n.cats, length(name_channel), length(index_state)),
    #                                      dimnames = list('params'=seq_len(n.cats), 'channel'=name_channel, 'state'=index_state)
    #                                    ) %>%
    #                                   apply(MARGIN=c(2,3), FUN=prob2diff, trunc=trunc) %>%
    #                                   array_tree(margin = c(3,2))
    #                                    }
    #                                  ),
    #

  #random starting values
    init.state.probs = map(n.states, ~dirichlet::rdirichlet(1, rep(1,.x)) %>% as.vector()), #random starting values
    #init.trans.P = map(n.states_cluster, ~map(.x, function(states){(dirichlet::rdirichlet(states, rep(1,states))+diag(states)) %>% {./rowSums(.)}}) %>% Matrix::bdiag() %>% as.matrix()), #random transition probability
    #init_trans.p = map(init_trans_probs, function(transmat){ (transmat + 0.01)/rowSums(transmat + 0.01)}), #random transition probability +small off-diagonal probs
    init.params.random_channel_state = pmap(list(index_state, name_channel, n.cats, trunc),
                                     function(index_state, name_channel, n.cats, trunc){
                                       dirichlet::rdirichlet(length(name_channel)*length(index_state), rep(1, n.cats)) %>%
                                       t() %>% unlist() %>% apply(2, FUN=prob2diff, trunc=trunc) %>%
                                       array(
                                         dim = c(n.cats-1, length(name_channel), length(index_state)),
                                         dimnames = list('params'=c("thres1","diff12","diff23","diff34","diff45","diff56"), 'channel'=name_channel, 'state'=index_state)
                                       ) %>% array_tree(margin = c(3,2))
                                       }
                                     ),


    hmodel.in = pmap(list(emissiondist, init.params.random_channel_state),
                  function(emissiondist, params_channel_state){
                    params_channel_state %>% map_depth(2, ~lift_dl(emissiondist)(.)) %>% map_depth(1, ~lift_dl(msm::hmmMV)(.))
                    }
                  ),

    hconstraint = pmap(list(index_cluster_state, n.channels,n.states),
                       function(index_cluster_state, n.channels, n.states){
                         hconstr <- rep((index_cluster_state-1)*n.channels, each=n.channels)+rep(seq_len(n.channels), n.states)
                         list(diff12=hconstr, diff23=hconstr, diff34=hconstr, diff45=hconstr, diff56=hconstr)
                      }
                  ),

    hranges = pmap(list(trunc, n.states, n.channels),
                   function(trunc, n.states, n.channels){
                     totalpars <- n.states*n.channels
                      list(thres1 = list(lower = rep(-trunc, totalpars), upper = rep(trunc, totalpars)),
                           diff12 = list(lower = rep(0, totalpars), upper = rep(trunc*2, totalpars)),
                           diff23 = list(lower = rep(0, totalpars), upper = rep(trunc*2, totalpars)),
                           diff34 = list(lower = rep(0, totalpars), upper = rep(trunc*2, totalpars)),
                           diff45 = list(lower = rep(0, totalpars), upper = rep(trunc*2, totalpars)),
                           diff56 = list(lower = rep(0, totalpars), upper = rep(trunc*2, totalpars))
                      )
                   }),

    n.pars = pmap_dbl(list(n.cats, n.channels, n.states_cluster, n.states, n.clusters),
                 function(n.cats, n.channels, n.states_cluster, n.states, n.clusters){
                   n.qpars <- sum(n.states_cluster*(n.states_cluster-1))
                   n.threspars <- n.states*n.channels
                   n.diffpars <- (n.cats-2)*n.clusters*n.channels
                   n.estinitprobs <- (n.states-1)
                   n.qpars+n.threspars+n.diffpars+n.estinitprobs
                   }
                )
)



# fit ordinal hmms from list of inputs in configuration ---------------------------------------------

## send this to iCSF as this will take a LONG time to run


library('furrr')
parallel::detectCores()
n_cores <- 10

clusters <- parallel::makeCluster(n_cores)
plan(cluster, workers=clusters)

tictoc::tic()
hmmslist_ord0 <-
  configuration_ord %>%
  mutate(
    hmmfit = future_pmap(
      list(chr.n.states_cluster, init.state.probs, init.trans.Q, hmodel.in, hconstraint, hranges, n.pars),
      function(chr.n.states_cluster, init.state.probs, init.trans.Q, hmodel.in, hconstraint, hranges, n.pars)
      {
        start_time <- Sys.time()
        mod<-msm::msm(
          formula = responseMV~fup_days,
          subject = ptid,
          data = data_wide,
          initprobs = init.state.probs,
          est.initprobs = TRUE,
          qmatrix = init.trans.Q,
          gen.inits = FALSE,
          hmodel = hmodel.in,
          hconstraint = hconstraint,
          hranges = hranges,
          obstype = 1,
          opt.method = "optim",
          control = list(fnscale=10000,
                         maxit=10000, #default=500 for Nelder mead
                         reltol=1e-10, #default~=1e-08
                         ndeps = rep(1e-05, n.pars) #default=1e-3
                         ),
           method="Nelder-Mead" #default="BFGS"
        )

      }#,.progress = TRUE
    ),
  )
tictoc::toc()
parallel::stopCluster(clusters)
plan(sequential)

write_rds(hmmslist_ord0, path = here::here("data", str_c("hmmslist0 ordinal grp seed is ", seed,".rds")), compress = "gz")

