
library('msm')


# process fitted hmm -------------------------------------------------------

## derive response probabilities for each state and combine lists

hmms_cat <- read_rds(path=here::here("data", "hmmslist0 categorical grp seed is 123.rds")) %>%
  mutate(
    distrtype = 'categorical',
    hmodellist = map(hmmfit, ~hmodel2list(.$hmodel, hmmdist=FALSE)),
    problist = map_depth(hmodellist, 3, ~unlist(.)[-8]),

    n.qpars = map_dbl(n.states_cluster, ~sum(.*(.-1))),
    n.initpars = (n.states-1),
    n.thetapars = n.states*n.channels*(n.cats-1)
  )


hmms_ord <- read_rds(path=here::here("data", "hmmslist0 ordinal grp seed is 123.rds")) %>%
  mutate(
    distrtype = 'ordinal',
    hmodellist = map(hmmfit, ~hmodel2list(.$hmodel, hmmdist=FALSE)),
    problist = map_depth(hmodellist, 3, ~diff2prob(unlist(.)[-7], 5)),

    n.qpars = map_dbl(n.states_cluster, ~sum(.*(.-1))),
    n.initpars = (n.states-1),
    n.thetapars = (n.states*n.channels) + ((n.cats-2)*n.clusters*n.channels)
  )


hmms0 <- bind_rows(hmms_cat, hmms_ord) %>%
  mutate(
    convergence = map_int(hmmfit,~.$opt$convergence),
    n.obs = map_dbl(hmmfit, ~ nrow(msm:::model.frame.msm(.))),
    #n.obs.size = map_dbl(hmmfit, ~ nrow(msm:::model.frame.msm(.))),
    n.responses = map_dbl(hmmfit, ~ sum(!is.na(select(msm:::model.frame.msm(.), starts_with("(state)"))$`(state)`)) ),

    n.estinitprobs = (n.states-1),


    logLik = map_dbl(hmmfit, logLik),
    AIC = map_dbl(hmmfit, AIC),
    AIC.manual = map2_dbl(logLik, n.pars, ~2*.y-2*.x),

    #true BIC is somewhere between these two values
    BIC = log(n.obs)*n.pars -(2*logLik),
    BIC.nonmiss = log(n.responses)*n.pars -(2*logLik),


    fittedQ = map(hmmfit, ~msm::qmatrix.msm(., ci="none")),
    fittedP1 = map(hmmfit, ~msm::pmatrix.msm(., ci = 'none', t=1)),
  )

hmms0 %>%
  select(distrtype, chr.n.states_cluster, n.states_cluster, n.states, n.clusters, n.obs, n.responses, n.qpars, n.initpars, n.thetapars, n.pars, AIC, BIC, BIC.nonmiss) %>%
  write_rds(path=here::here("data", "model_performance.rds"))


## compare BIC of moedls
plot_BIC <-
  hmms0 %>%
  select(chr.n.states_cluster, n.states, n.clusters, distrtype, AIC, BIC, BIC.nonmiss) %>%
  #spread(distrtype, BIC) %>%
  ggplot() +
  geom_abline(aes(intercept=0, slope=1), linetype='dashed')+
  geom_point(aes(x=BIC, y=BIC.nonmiss, colour=as.factor(distrtype))) +
  ggrepel::geom_text_repel(aes(x=BIC, y=BIC.nonmiss, label = str_replace_all(chr.n.states_cluster, "_", "")))+
  coord_fixed()+
  theme_bw()+
  NULL
ggsave(filename=here::here("figures",str_c("BIC ",".png")), plot_BIC, units='cm', height=30, width=40)


## compare BIC of moedls
plot_BIC <-
  hmms0 %>%
  select(chr.n.states_cluster, distrtype, BIC) %>%
  spread(distrtype, BIC) %>%
  ggplot() +
  geom_abline(aes(intercept=0, slope=1), linetype='dashed')+
  geom_point(aes(x=categorical, y=ordinal))+#, colour=as.factor(distrtype))) +
  ggrepel::geom_text_repel(aes(x=categorical, y=ordinal, label = str_replace_all(chr.n.states_cluster, "_", "")))+
  coord_fixed()+
  theme_bw()+
  NULL





####

# model ordinal 2_2_3 has the lowest BIC, so we use this model for our simulations.



hmms0 %>%
  filter(
    chr.n.states_cluster=="2_2_3",
    distrtype == "ordinal"
  ) %>%
  write_rds(path = here::here("data", "hmm_ord_2_2_3.rds"), compress="gz")





hmm_single <- read_rds(path = here::here("data", "hmm_ord_2_2_3.rds"))

hmms_evaluate <- hmm_single  %>%
  mutate(


    states_subject_time =
      pmap(list(hmmfit, index_state, index_cluster_state),
          function(hmmfit, index_state, index_cluster_state){
            msm::viterbi.msm(hmmfit, normboot=FALSE) %>%
            mutate(
              index_state = list(index_state),
              index_cluster_state = list(index_cluster_state),
              state_prob = array_tree(.$pstate, 1),
              #n.nans = map_dbl(state_prob, ~sum(is.nan(.))),
              state_prob = map(state_prob, function(sp){spnew<-sp; spnew[is.nan(sp)]<-1-sum(sp, na.rm=TRUE); spnew}),#lot's of NaNs - assume other values are correct and impute
              cluster_prob = map2(state_prob, index_cluster_state, function(probs, cluster){aggregate(probs, list(a=cluster), FUN=sum)$x}),
              guess_state = fitted
              ) %>%
            select(-fitted, -pstate)
          }
      ),

    pe_probs_state =
      pmap(
        list(n.states, n.cats, index_state, index_cluster_state, index_channel, name_channel, problist),
        function(n.states, n.cats, index_state, index_cluster_state, index_channel, name_channel, problist)
        {
          tibble(index_state, index_cluster_state, problist) %>%
          unnest() %>%
          mutate(
            name_channel = rep(name_channel, n.states),
            index_channel = rep(index_channel, n.states),
            response = list(seq_len(n.cats))
          ) %>%
          rename(state=index_state, cluster=index_cluster_state, probs=problist) %>%
          arrange(state, cluster, index_channel) %>%
          mutate(
            cml.probs=map(probs, cumsum),
            #sum.probs=map_dbl(probs, sum) #UNIT TEST should be 1
          ) %>%
          select(state, cluster, index_channel, name_channel, response, probs, cml.probs)
      }
    ),


    pe_obs_state =
      pmap(
        list(states_subject_time, name_channel, index_channel),
        function(states_subject_time, name_channel, index_channel)
        {
          states_subject_time %>%
          select(-cluster_prob) %>%
          unnest() %>%
          group_by(index_state, index_cluster_state) %>%
          mutate(
            state_prob_sum=sum(state_prob),
            #response0 = list(seq_len(n.cats))
          ) %>%
          gather("item", "response", starts_with("observed.")) %>%
          rename(state=index_state, cluster=index_cluster_state) %>%
          group_by(state, cluster, item) %>%
          mutate(
            name_channel=str_remove(item, "observed."),
            #state_weight2=sum(state_prob),
            response0=factor(ifelse(is.na(response),0,response), levels=c(0, seq_len(7))) #convert NAs to zero and create factor so that grouping doesn't drop zero frequencies
          ) %>%
          group_by(state, cluster, name_channel, response0, .drop=FALSE) %>%
          summarise(
            obs_freq0=sum(if_else(state_prob==0, 0, state_prob/state_prob_sum)) # this is the frequency of responses, weighted by the forward-backward probability of being in each state
          ) %>%
          group_by(state, cluster, name_channel) %>%
          summarise(
            response0=list(response0),
            obs_freq0=list(obs_freq0)
          ) %>%
          mutate(
            obs_freq=map(obs_freq0, ~.[-1]/sum(.[-1])),
            cml.obs_freq=map(obs_freq, cumsum),
            #freqsum0=map_dbl(obs_freq0, sum), # UNIT TEST should be 1
            #freqsum=map_dbl(obs_freq, sum) # UNIT TEST should 1
          ) %>%
          left_join(tibble(name_channel, index_channel), by=c('name_channel')) %>%
          select(state, cluster, index_channel, name_channel, response0, obs_freq0, obs_freq, cml.obs_freq)
        }
      ),

    state_weights =
      pmap(list(states_subject_time),
           function(states_subject_time)
           {
             states_subject_time %>%
             select(-cluster_prob) %>%
             unnest() %>%
             rename(state=index_state, cluster=index_cluster_state) %>%
             group_by(state, cluster, subject) %>%
             summarise(
               state_prob_subject=mean(state_prob),
               #state_prob_subject2 = log(prod(exp(state_prob)))/n()
               ) %>%
             ungroup()
           }
      ),

     pe_obs_state_subject =
      pmap(
        list(states_subject_time, state_weights, name_channel, index_channel),
        function(states_subject_time, state_weights, name_channel, index_channel)
        {
          states_subject_time %>%
            select(-cluster_prob) %>%
            unnest() %>%
            group_by(index_state, index_cluster_state, subject) %>%
            mutate(
              state_prob_sum=sum(state_prob)
            ) %>%
            gather("item", "response", starts_with("observed.")) %>%
            rename(state=index_state, cluster=index_cluster_state) %>%
            group_by(state, cluster, item) %>%
            mutate(
              name_channel=str_remove(item, "observed."),
              state_weight=sum(state_prob),
              response0=factor(ifelse(is.na(response),0,response), levels=c(0, seq_len(7))) #convert NAs to zero
            ) %>%
            group_by(state, cluster, subject, name_channel, response0, .drop=FALSE) %>%
            summarise(
              obs_freq0=sum(if_else(state_prob==0, 0, state_prob/state_prob_sum)),# this is the frequency of responses, weighted by the forward probability of being in each state
            ) %>%
            group_by(state, cluster, subject, name_channel) %>%
            summarise(
              response0=list(response0),
              obs_freq0=list(obs_freq0),
            ) %>%
            mutate(
              obs_freq=map(obs_freq0, ~.[-1]/sum(.[-1])),
              cml.obs_freq=map(obs_freq, cumsum),
              #freqsum0=map_dbl(obs_freq0, sum), # UNIT TEST should be 1
              #freqsum=map_dbl(obs_freq, sum) # UNIT TEST should 1
            ) %>%
            left_join(tibble(name_channel, index_channel), by=c('name_channel')) %>%
            left_join(state_weights) %>%
            select(state, cluster, ptid=subject, index_channel, name_channel, response0, obs_freq0, obs_freq, cml.obs_freq, state_prob_subject)
        }
      ),

    pe_eval_state=
      pmap(
        list(pe_probs_state, pe_obs_state),
        ~left_join(.x, .y, by=c('state', 'cluster', 'index_channel', 'name_channel')) %>%
          mutate(
            wasserstein=
              pmap_dbl(
                list(response, obs_freq, probs),
                    function(response, obs_freq, probs){wassdiscrete(obs_freq, probs, norm=1, support=response)}
              ),
            wdfit = 1-(wasserstein/(map_int(response, length)-1))
          )
      ),

    pe_eval_state_subject=
      pmap(
        list(pe_probs_state, pe_obs_state_subject),
        ~left_join(.x, .y, by=c('state', 'cluster', 'index_channel', 'name_channel')) %>%
          mutate(

            wasserstein =
                pmap_dbl(
                  list(response, obs_freq, probs),
                      function(response, obs_freq, probs){wassdiscrete(obs_freq, probs, norm=1, support=response)}
                ),
            wdfit = 1-(wasserstein/(map_int(response, length)-1))
          )
        ),

  )


## plots



library('Cairo')
setHook(packageEvent("grDevices", "onLoad"),
        function(...) grDevices::X11.options(type="cairo")
)
options(device="x11")
options(bitmapType='cairo')


## pp plots
plot_ppstep <-
  hmms_evaluate %>%
  #filter(chr.n.states_cluster==config) %>%
  `[[`('pe_eval_state') %>% `[[`(1) %>%
  mutate(
    cml.probs=map(cml.probs,~c(0,.)),
    cml.obs_freq=map(cml.obs_freq,~c(0,.))
  ) %>%
  unnest(response0, cml.probs, cml.obs_freq) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1), colour='grey')+
  geom_step(aes(x=cml.probs, y=cml.obs_freq), size=0.1, colour='red') +
  coord_fixed()+
  scale_x_continuous(name="cumulative expected response", limits=c(0,1))+
  scale_y_continuous(name="cumulative observed response", limits=c(0,1))+
  lemon::facet_rep_grid(name_channel~str_c(LETTERS[cluster],"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("ppstep ", distr," ",config,".png")), plot_ppstep, units='cm', height=20, width=30)




plot_ppline <-
  hmms_evaluate %>%
  #filter(chr.n.states_cluster==config) %>%
  `[[`('pe_eval_state_subject') %>% `[[`(1) %>%
  mutate(
    cml.probs=map(cml.probs,~c(0,.)),
    cml.obs_freq=map(cml.obs_freq,~c(0,.))
  ) %>%
  unnest(response0, cml.probs, cml.obs_freq) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1), colour='grey')+
  geom_line(aes(x=cml.probs, y=cml.obs_freq, group=ptid, alpha=state_prob_subject), size=0.1, colour='red') +
  coord_fixed()+
  scale_x_continuous(name="cumulative expected response", limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_y_continuous(name="cumulative observed response", limits=c(0,1), breaks=c(0, 0.5, 1))+
  scale_alpha_continuous(range = c(0, 1))+
  lemon::facet_rep_grid(name_channel~str_c(LETTERS[cluster],"\n",state))+
  theme_bw(base_size=13)+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("ppline pp ", distr," ",config,".png")), plot_ppline, units='cm', height=15, width=30)




## qq plots


plot_response_qqline <-
  hmms_evaluate %>%
  `[[`('pe_eval_state_subject') %>% `[[`(1) %>%
  filter(state_prob_subject>0.0000001) %>%
  mutate(
    p_in=list(seq(0,1,0.01)),
    q_obs=map2(obs_freq, p_in, ~qmultinom_itpl(.y, .x, is.prob.cumul = FALSE)),
    q_exp=map2(probs, p_in, ~qmultinom_itpl(.y, .x, is.prob.cumul = FALSE))
  ) %>%
  unnest(p_in, q_obs, q_exp) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1), colour='grey')+
  geom_step(aes(x=q_exp, y=q_obs, group=ptid, alpha=state_prob_subject), size=0.1, colour='red') +
  #geom_jitter(aes(x=q_exp, y=q_obs, group=ptid), alpha=0.1, size=0.1, width=0.3, height=0.3, colour='red') +
  coord_fixed()+
  scale_x_continuous(name="expected response", limits=c(0,7.5), breaks=1:7)+
  scale_y_continuous(name="observed response", limits=c(0,7.5), breaks=1:7)+
  lemon::facet_rep_grid(name_channel~str_c(cluster,"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("qqline pp ", distr," ",config,".png")), plot_response_qqline, units='cm', height=20, width=30)


plot_response_qqpoint <-
  hmms_evaluate%>%
  `[[`('pe_eval_state_subject') %>% `[[`(1) %>%
  mutate(
    p_in=list(seq(0,1,0.01)),
    q_obs=map2(obs_freq, p_in, ~qmultinom(.y, unname(.x), is.prob.cumul = FALSE)),
    q_exp=map2(probs, p_in, ~qmultinom(.y, unname(.x), is.prob.cumul = FALSE))
  ) %>%
  unnest(p_in, q_obs, q_exp) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1), colour='grey')+
  geom_jitter(aes(x=q_exp, y=q_obs, group=ptid,alpha=state_prob_subject), size=0.1, width=0.3, height=0.3, colour='red') +
  coord_fixed()+
  scale_x_continuous(name="expected response", limits=c(0,7.5), breaks=1:7)+
  scale_y_continuous(name="observed response", limits=c(0,7.5), breaks=1:7)+
  lemon::facet_rep_grid(name_channel~str_c(cluster,"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("qqpoint pp ", distr," ",config,".png")), plot_response_qqpoint, units='cm', height=20, width=30)




## density plots

plot_response_probs <-
  hmms_evaluate %>%
  `[[`('pe_eval_state') %>% `[[`(1) %>%
  unnest(response, probs) %>%
  ggplot()+
  geom_bar(aes(x=factor(response), y=probs, fill=factor(response)), stat='identity') +
  scale_x_discrete(name="expected response", breaks=1:7, labels=c("1":"7"), drop=FALSE)+
  scale_y_continuous(name="frequency", expand=c(0,0), limits=c(0,1))+
  scale_fill_manual(values=c(viridisLite::viridis(n=7)))+
  lemon::facet_rep_grid(name_channel~str_c(LETTERS[cluster],"\n",state))+
  theme_bw(base_size=13)+
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    axis.text.x=element_blank(), axis.ticks = element_blank(),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.border = element_blank(),
    axis.line.x=element_line(colour='black'),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("dens prob ", distr," ",config,".png")), plot_response_probs, units='cm', height=15, width=25)


plot_response_obs <-
  hmms_evaluate %>%
  `[[`('pe_eval_state') %>% `[[`(1) %>%
  #select(-obs_freq0) %>%
  unnest(response0, obs_freq0) %>%
  ggplot()+
  geom_bar(aes(x=factor(response0), y=obs_freq0,  fill=factor(response0)), stat='identity') +
  scale_x_discrete(name="observed response", breaks=1:7, labels=c("1":"7"), drop=FALSE)+
  scale_y_continuous(name="frequency", expand=c(0,0), limits=c(0,1))+
  scale_fill_manual(values=c("grey",viridisLite::viridis(n=7)))+
  lemon::facet_rep_grid(name_channel~str_c(LETTERS[cluster],"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    axis.text.x=element_blank(), axis.ticks = element_blank(),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.border = element_blank(),
    axis.line.x=element_line(colour='black'),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("dens obs ", distr," ",config,".png")), plot_response_obs, units='cm', height=20, width=30)


plot_wd <-
  hmms_evaluate %>%
  `[[`('pe_eval_state_subject') %>% `[[`(1) %>%
  filter(state_prob_subject>0.0000001) %>%
  ggplot()+
  #geom_jitter(aes(x=state, y=wdfit), width=0.2, alpha=0.2) +
  # geom_violin(aes(x=state, y=wdfit), width=0.2, alpha=0.2) +
  geom_dotplot(aes(x=1, y=wdfit, alpha=state_prob_subject), binaxis='y', binwidth=0.02, stackdir='center', width=0.2) +
  scale_x_continuous(name=" ")+
  scale_y_continuous(name="1 - normalised Wasserstein distance", limits=c(0,1))+
  lemon::facet_rep_grid(name_channel~str_c(cluster,"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right = element_text(angle=0, vjust=0.05, hjust=1),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    strip.text.y = element_text(angle=0, hjust=0.5),
    #panel.border = element_blank(),
    axis.line.x=element_line(colour='black'),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL
ggsave(filename=here::here("figures",str_c("wd ", distr," ",config,".png")), plot_wd, units='cm', height=20, width=30)









compare_wd <-
  left_join(
    hmm_evaluationcat %>%
    select(chr.n.states_cluster, n.states, n.clusters, pe_eval_state_subject) %>%
    mutate(
      pe_eval_state_subject = map(pe_eval_state_subject, . %>% select(state, cluster, index_channel, name_channel, ptid, wasserstein, wdfit))
      ) %>%
    unnest() %>%
    rename(wasserstein.cat = wasserstein, wdfit.cat = wdfit),
    hmm_evaluationord %>%
      select(chr.n.states_cluster, n.states, n.clusters, pe_eval_state_subject) %>%
      mutate(
        pe_eval_state_subject = map(pe_eval_state_subject, . %>% select(state, cluster, index_channel, name_channel, ptid, wasserstein, wdfit))
      ) %>%
      unnest() %>%
      rename(wasserstein.ord = wasserstein, wdfit.ord = wdfit)
  )


compare_wd %>%
  filter(chr.n.states_cluster=="2_2") %>%
  ggplot() +
  geom_abline(aes(intercept=0, slope=1), linetype='dashed')+
  geom_point(aes(x=wdfit.cat, y=wdfit.ord), binaxis='y', binwidth=0.02, stackdir='center', width=0.2) +
  coord_fixed()+
  #scale_x_continuous(name=" ")+
  #scale_y_continuous(name="1 - normalised Wasserstein distance", limits=c(0,1))+
  lemon::facet_rep_grid(name_channel~str_c(cluster,"\n",state))+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right = element_text(angle=0, vjust=0.05, hjust=1),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    strip.text.y = element_text(angle=0, hjust=0.5),
    #panel.border = element_blank(),
    axis.line.x=element_line(colour='black'),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL


























View(data_CT_encounter_itemgrp)
"P01002" "P01004" "P01005" "P01007" "P01011" "P01024" "P01027" "P01036" "P01041" "P01043" "P01044" "P01045" "P01054" "P01066" "P01068"
[16] "P01072" "P01079" "P01081" "P02010" "P02014" "P02015" "P02018" "P02021" "P02025" "P02026" "P02029" "P02031" "P02037" "P02040" "P02046"
[31] "P02050" "P02052" "P02055" "P02058" "P02059" "P02061" "P02062" "P02065" "P02071" "P02077" "P10"    "P103"   "P104"   "P110"   "P112"
[46] "P116"   "P117"   "P119"   "P120"   "P127"   "P131"   "P133"

ptid_exp <- "P02055"
example_traj <- data_CT_encounter_itemgrp %>%
  filter(ptid %in% ptid_exp) %>%
  #filter(ptid %ni% c("P02029","P02065", "P104")) %>% #,
  ggplot() +
  facet_grid(rows=vars(itemgrp))+
  geom_point(aes(x=fup_days, y=score, colour=as.factor(score)), size=0.8) +
  labs(colour="", x="days", y="response")+
  scale_colour_manual(values=c(viridisLite::viridis(n=7)), na.value="transparent", breaks=1:7, labels=c(1:7))+
  scale_x_continuous(limits=c(0,90), breaks=seq(0,100,7))+
  scale_y_continuous(breaks=c(1,4,7), limits=c(0,8))+
  labs(x="Day", y=NULL, fill="Response")+
  theme_bw()+
  theme(
    legend.position="bottom",legend.direction="horizontal",
    strip.background = element_blank(),
    strip.text.y = element_text(angle=0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  guides(colour=guide_legend(nrow=1))

ggsave(filename=here::here("figures",str_c("traj ",ptid_exp,".png")), example_traj, units='cm', height=10, width=20)

data_CT_encounter_itemgrp %>%
  filter(ptid %in% c("P01004")) %>%
  #left_join(pe_cluster_probs %>% filter(n_states_chr=="1_2_3_3") %>% group_by(ptid) %>% summarise(cluster_guess=first(cluster_guess))) %>%
  #filter(ptid %ni% c("P02029","P02065", "P104")) %>% #,
  ggplot() +
  #facet_grid(rows=vars(itemgrp))+
  geom_tile(aes(x=fup_days, y=fct_rev(itemgrp), fill=as.factor(score)), height=0.6, width=0.1, colour='transparent') +
  labs(fill="", x="days", y="response")+
  scale_fill_manual(values=c(viridisLite::viridis(n=7)), na.value="transparent", breaks=1:7, labels=c(1:7))+
  scale_x_continuous(limits=c(0,100))+
  labs(x="Day", y=NULL, fill="Response")+
  theme_bw()+
  theme(
    legend.position="bottom",legend.direction="horizontal",
    strip.text = element_blank(),
    strip.text.y = element_text(angle=0),
    panel.grid.major.y = element_blank()
  )+
  guides(fill=guide_legend(nrow=1))






plotday <-
  data_CT_encounter_itemgrp %>%
  #left_join(pe_cluster_probs %>% filter(n_states_chr=="1_2_3_3") %>% group_by(ptid) %>% summarise(cluster_guess=first(cluster_guess))) %>%
  #filter(ptid %ni% c("P02029","P02065", "P104")) %>% #,
  ggplot() +
  geom_tile(aes(x=fup_days, y=itemgrp, fill=as.factor(score)), height=0.6, width=1, colour="transparent") +
  #geom_text(aes(x=89, y="Anxiety", label=LETTERS[cluster_guess]), hjust=0, colour='red')+
  facet_wrap(~ptid, ncol=10)+
  scale_fill_manual(values=c(viridisLite::viridis(n=7)), na.value="transparent", breaks=1:7)+
  scale_x_continuous(limits=c(0,100))+
  labs(x="Day", y=NULL, fill="Response")+
  theme_bw()+
  theme(
    legend.position="bottom",legend.direction="horizontal",
    strip.text = element_blank(),
    panel.grid.major.y = element_blank()
  )+
  guides(fill=guide_legend(nrow=1))
ggsave(filename="figures/mean score per day.png", fig_meanscore_day, units='cm', height=20, width=30)
