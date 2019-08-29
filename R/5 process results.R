#library('tidyverse')
#library('msm')
#import::from(magrittr, "%$%")
#import::from(willsutils, "%ni%")


data_canon <- read_rds(here::here("data", "data_canon.rds")) %>%
  mutate(itemgrp = factor(itemgrp, levels=c("Anxiety", "Delusions", "Depression", "Hallucinations")))

data_wide <- read_rds(here::here("data", "data_wide.rds"))


# response distribution ---------------------------------------------------

hmm_223 <- read_rds(path = here::here("data", "hmm_ord_2_2_3.rds"))

hmm_obj <- hmm_223 %>%
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

    probstbl = pmap(
      list(index_state, index_cluster_state, name_channel, problist),
      function(index_state, index_cluster_state, name_channel, problist){
        dat <-

          tibble(probs = problist) %>%
          mutate(
            state=index_state,
            cluster=index_cluster_state,
            statename=str_c("C",cluster,"S",ave(parse_number(state), cluster, FUN=function(x){x-min(x)+1}))
          ) %>%
          unnest() %>%
          mutate(
            item = rep(name_channel, length(index_state)),
            response = map(probs, ~seq_along(.))
          ) %>%
          unnest()
      }
    )
  )


# plot chosen model -------------------------------------------------------


plot_respdistr <-
  pluck(hmm_obj, "probstbl", 1) %>%
  ggplot()+
  geom_bar(aes(x=factor(response), y=probs, fill=factor(response)), stat='identity') +
  scale_x_discrete(name="Response", breaks=c(1,4,7), labels=c("1","4","7"), drop=FALSE)+
  scale_y_continuous(name="Density", expand=c(0,0),breaks=c(0,0.5,1), limits=c(0,1))+
  scale_fill_manual(values=c(viridisLite::viridis(n=7)))+
  lemon::facet_rep_grid(item ~ statename)+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    #strip.text = element_blank(),
    strip.background = element_blank(),
    axis.title.y.right= element_text(angle = 0, vjust=0.05, hjust=1),
    #axis.text.x=element_blank(),
    axis.ticks = element_blank(),
    strip.text.y = element_text(angle=0, hjust=0.5),
    panel.border = element_blank(),
    axis.line.x=element_line(colour='black'),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )+
  NULL

ggsave(here::here("ord_2_2_3_results", "fig_respdistr.png"), plot=plot_respdistr)




# plot PP plots -------------------------------------------------------



hmm_evaluation <-
  hmm_obj %>%
  mutate(

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


    pe_eval_state_subject=
      pmap(
        list(pe_probs_state, pe_obs_state_subject),
        ~left_join(.x, .y,  by=c('state', 'cluster', 'index_channel', 'name_channel')) %>%
          mutate(
            statename=str_c("C",cluster,"S",ave(parse_number(state), cluster, FUN=function(x){x-min(x)+1})),
            wasserstein =
              pmap_dbl(
                list(response, obs_freq, probs),
                function(response, obs_freq, probs){wassdiscrete(obs_freq, probs, norm=1, support=response)}
              ),
            wdfit = 1-(wasserstein/(map_int(response, length)-1))
          )
      ),
  )

plot_ppline <-
  pluck(hmm_evaluation, "pe_eval_state_subject", 1) %>%
  mutate(
    cml.probs=map(cml.probs,~c(0,.)),
    cml.obs_freq=map(cml.obs_freq,~c(0,.))
  ) %>%
  unnest(response0, cml.probs, cml.obs_freq) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1), colour='grey', linetype='dotted')+
  geom_line(aes(x=cml.probs, y=cml.obs_freq, group=ptid, alpha=state_prob_subject), size=0.1, colour='#443A83FF') +
  coord_fixed()+
  scale_x_continuous(name="cumulative expected response", limits=c(0,1), breaks=c(0, 0.5, 1), labels=c("0", "1/2", "1"))+
  scale_y_continuous(name="cumulative observed response", limits=c(0,1), breaks=c(0, 0.5, 1), labels=c("0", "1/2", "1"))+
  scale_alpha_continuous(range = c(0, 0.9))+
  lemon::facet_rep_grid(name_channel~statename)+
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
ggsave(here::here("ord_2_2_3_results", "fig_ppline.png"), plot=plot_ppline)




plot_qqline <-
  pluck(hmm_evaluation, "pe_eval_state_subject", 1) %>%
  filter(state_prob_subject>0.0000001) %>%
  mutate(
    p_in=list(seq(0,1,0.01)),
    q_obs=map2(obs_freq, p_in, ~qmultinom_itpl(.y, .x, is.prob.cumul = FALSE)),
    q_exp=map2(probs, p_in, ~qmultinom_itpl(.y, .x, is.prob.cumul = FALSE))
  ) %>%
  unnest(p_in, q_obs, q_exp) %>%
  ggplot()+
  geom_abline(aes(intercept=0, slope=1),  colour='grey', linetype='dotted')+
  geom_step(aes(x=q_exp, y=q_obs, group=ptid, alpha=state_prob_subject), size=0.1, colour='#443A83FF') +
  #geom_jitter(aes(x=q_exp, y=q_obs, group=ptid), alpha=0.1, size=0.1, width=0.3, height=0.3, colour='red') +
  coord_fixed()+
  scale_x_continuous(name="expected response", limits=c(0,7.5), breaks=1:7)+
  scale_y_continuous(name="observed response", limits=c(0,7.5), breaks=1:7)+
  scale_alpha_continuous(range = c(0, 1))+
  lemon::facet_rep_grid(name_channel~statename)+
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
ggsave(here::here("ord_2_2_3_results", "fig_qqline.png"), plot=plot_qqline)



# simulation results ------------------------------------------------------


results1 <- readRDS(file=here::here("data", "simulations", "simresults_223_ss1_bs2.rds"))
results2 <- readRDS(file=here::here("data", "simulations", "simresults_223_ss3_bs4.rds"))
results3 <- readRDS(file=here::here("data", "simulations", "simresults_223_ss5_bs6.rds"))

results <- bind_rows(
  results1,
  results2,
  results3
  ) %>%
  mutate(
    cluster=case_when(startstate==1 ~ "C1",
                      startstate==3 ~ "C2",
                      startstate==5 ~ "C3",
                      )
  )

findclosest <-function(needed, avail){

  tempfunc <- Vectorize(function(needed, avail){
    if( !identical(sort(avail), avail) ) stop("avail not in ascending order")

   # if( length(needed) == 0 | is.na(needed))
  #    NA_real_
  #  else
      avail[which(avail>=needed)[1]]

  }, vectorize.args="needed")

  x <- tempfunc(needed, avail)
  if( length(x)==0) NA_real_
  else x
}



evaluateavg <- results %>%
  #group_by(tau, repID, maxtime) %>%
  transmute(
    tau, startstate, badstate, repID, maxtime, cluster,

    n=map_int(adp, nrow),
    #badstate_ttime,
    #badstate_ttime2 = map(badstate_ttime, ~function(x) {if length(.)==0, NA_real_, .)),

    meandelay = map2_dbl(adp, badstate_ttime, ~mean(findclosest(.y, .x$time_J)-.y, na.rm=TRUE)),
    #mindelay = map2_dbl(adp, badstate_ttime, ~min(findclosest(.y, .x$time_J)-.y, na.rm=TRUE)),
    #maxdelay = map2_dbl(adp, badstate_ttime, ~max(findclosest(.y, .x$time_J)-.y, na.rm=TRUE)),
    meaninterval = map_dbl(adp, ~mean(diff(.x$time_J), na.rm=TRUE)),

    obsfreq = pmap_dbl(list(maxtime, adp),
                        function(maxtime, adp){
                          nrow(adp)/maxtime
                        }
                        ),


    obsfreqalert = pmap_dbl(list(state_ttime, state_chain, badstate, maxtime, adp),
                        function(state_ttime, state_chain, badstate, maxtime, adp){
                          timeinstate <- sum(diff(c(state_ttime, maxtime))[state_chain==badstate])
                          obsinstate <- sum(state_chain[findInterval(adp$time_J, state_ttime)]==badstate)
                          min(obsinstate/timeinstate,4)
                        }
                        ),

    obsfreqnotalert = pmap_dbl(list(state_ttime, state_chain, badstate, maxtime, adp),
                        function(state_ttime, state_chain, badstate, maxtime, adp){
                          timeinstate <- sum(diff(c(state_ttime, maxtime))[state_chain!=badstate])
                          obsinstate <- sum(state_chain[findInterval(adp$time_J, state_ttime)]!=badstate)
                          min(obsinstate/timeinstate,4)
                        }
                        ),
  )

evaluate <- results %>%
  #group_by(tau, repID, maxtime) %>%
  transmute(
    tau, startstate, badstate, repID, maxtime, cluster,

    badstatetrans = map_int(badstate_ttime, ~length(.x)),
    n=map_int(adp, nrow),
    delay = map2(adp, badstate_ttime, ~findclosest(.y, .x$time_J)-.y),
  ) %>%
  unnest()





plot_freq <-
  evaluateavg %>%
  filter(tau!=0) %>%
  #filter(cluster=="C1") %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             #cols=vars(tau),
             scales="free_y"
             )+
  geom_density(aes(x=obsfreq, colour=as.factor(tau)), size=1, alpha=1)+
  #geom_histogram(aes(x=n/100, fill=as.factor(tau)), size=1, alpha=0.3)+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 4.1)) +
  #scale_color_brewer(type='seq', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6)[-1], na.value="transparent", labels=c("2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Mean observation frequency (obs/day)', colour='tau')+
  theme_bw()+
  theme(legend.position='none')

ggsave(here::here("ord_2_2_3_results", "sim_obsfreq.png"), plot=plot_freq, width=20, height=15, units="cm")


plot_freqalert <-
  evaluateavg %>%
  filter(tau!=0) %>%
  #filter(cluster=="C1") %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             #cols=vars(tau),
             scales="free_y"
             )+
  geom_density(aes(x=obsfreqalert, colour=as.factor(tau)), size=1, alpha=1)+
  #geom_histogram(aes(x=n/100, fill=as.factor(tau)), size=1, alpha=0.3)+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 4.1)) +
  #scale_color_brewer(type='seq', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6)[-1], na.value="transparent", labels=c("2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Mean observation frequency in alert state (obs/day)', colour='tau')+
  theme_bw()+
  theme(legend.position='none')

ggsave(here::here("ord_2_2_3_results", "sim_obsfreqalert.png"), plot=plot_freq, width=20, height=15, units="cm")


plot_freqnotalert <-
  evaluateavg %>%
  filter(tau!=0) %>%
  #filter(cluster=="C1") %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             #cols=vars(tau),
             scales="free_y"
             )+
  geom_density(aes(x=obsfreqnotalert, colour=as.factor(tau)), size=1, alpha=1)+
  #geom_histogram(aes(x=n/100, fill=as.factor(tau)), size=1, alpha=0.3)+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 4.1)) +
  #scale_color_brewer(labels=c("2^-5", "2^-4","2^-3","2^-2","2^-1"), type='seq', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6)[-1], na.value="transparent", labels=c("2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Mean observation frequency not in alert state (obs/day)', colour='tau')+
  theme_bw()+
  theme(legend.position='right')

ggsave(here::here("ord_2_2_3_results", "sim_obsfreqnotalert.png"), plot=plot_freq, width=20, height=15, units="cm")


freqplots <- cowplot::plot_grid(plot_freq, plot_freqalert, plot_freqnotalert, ncol=1, axis="rl", align="hv")

ggsave(here::here("ord_2_2_3_results", "freqplots.png"), plot=freqplots, width=20, height=20, units="cm")



plot_avgdelay <-
  evaluateavg %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             scales="free_y")+
  geom_density(aes(x=meandelay, colour=as.factor(tau)), size=1, alpha=1)+
  scale_y_continuous(expand = c(0, 0))+
  #scale_color_brewer(labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), type='seq', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Mean delay before observation in alert state (days)', colour='tau')+
  theme_bw()

ggsave(here::here("ord_2_2_3_results", "sim_avgdelay.png"), plot=plot_avgdelay, width=20, height=15, units="cm")










plot_delay <-
  evaluate %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             scales='free_y')+
  geom_density(aes(x=delay, colour=as.factor(tau)), size=1, alpha=1)+
  scale_y_continuous(expand = c(0, 0)) +
  #scale_color_brewer(type='seq', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Delay (days)', colour='tau')+
  theme_bw()

ggsave(here::here("ord_2_2_3_results", "sim_delay.png"), plot=plot_delay, width=20, height=15, units="cm")


plot_interval <-
  evaluateavg %>%
  #filter(tau!=0) %>%
  #filter(cluster=="C1") %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             scales='free_y')+
  geom_density(aes(x=meaninterval, colour=as.factor(tau)), size=1)+
  scale_y_continuous(expand = c(0, 0)) +
  #scale_color_brewer(type='qual', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Mean interval between observation times', colour='tau')+
  theme_bw()

ggsave(here::here("ord_2_2_3_results", "sim_interval.png"), plot=plot_interval, width=20, height=15, units="cm")



plot_nbadtransitions <-
  evaluate %>%
  ggplot() +
  facet_grid(rows=vars(cluster),
             scales="free_y")+
  geom_bar(aes(x=badstatetrans, fill=as.factor(tau)), size=1, position ='dodge')+
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_brewer(type='qual', direction=-1)+
  scale_fill_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x='Transitions to alert state', colour='tau')+
  theme_bw()

ggsave(here::here("ord_2_2_3_results", "sim_transitions.png"), plot=plot_nbadtransitions, width=20, height=15, units="cm")




results1table <- results %>%
  #group_by(tau, repID, maxtime) %>%
  #filter(cluster=="C1") %>%
  transmute(
    tau, startstate, badstate, repID, maxtime, cluster,

    n=map_int(adp, nrow),

    #badstate_ttime,
    #badstate_ttime2 = map(badstate_ttime, ~function(x) {if length(.)==0, NA_real_, .)),
    badstatetrans = map_int(badstate_ttime, ~length(.x)),
    interval = map_dbl(adp, ~mean(diff(.$time_J), na.rm=TRUE)),
    delay = map2(adp, badstate_ttime, ~findclosest(.y, .x$time_J)-.y),


    obsfreqalert = pmap_dbl(list(state_ttime, state_chain, badstate, maxtime, adp),
                        function(state_ttime, state_chain, badstate, maxtime, adp){
                          timeinstate <- sum(diff(c(state_ttime, maxtime))[state_chain==badstate])
                          obsinstate <- sum(state_chain[findInterval(adp$time_J, state_ttime)]==badstate)
                          obsinstate/timeinstate
                        }
                        ),

    obsfreqnotalert = pmap_dbl(list(state_ttime, state_chain, badstate, maxtime, adp),
                        function(state_ttime, state_chain, badstate, maxtime, adp){
                          timeinstate <- sum(diff(c(state_ttime, maxtime))[state_chain!=badstate])
                          obsinstate <- sum(state_chain[findInterval(adp$time_J, state_ttime)]!=badstate)
                          obsinstate/timeinstate
                        }
                        ),

  ) %>%
  unnest() %>%
  group_by(tau, cluster) %>%
  summarise(

    n_Q1 = quantile(n, 0.25),
    n_Q2 = quantile(n, 0.5),
    n_Q3 = quantile(n, 0.75),

    frequency_Q1 = quantile(n/100, 0.25),
    frequency_Q2 = quantile(n/100, 0.5),
    frequency_Q3 = quantile(n/100, 0.75),

    interval_Q1 = quantile(interval, 0.25, na.rm=TRUE),
    interval_Q2 = quantile(interval, 0.5, na.rm=TRUE),
    interval_Q3 = quantile(interval, 0.75, na.rm=TRUE),

    delay_Q1 = quantile(delay, 0.25, na.rm=TRUE),
    delay_Q2 = quantile(delay, 0.5, na.rm=TRUE),
    delay_Q3 = quantile(delay, 0.75, na.rm=TRUE),

    freqalert_Q1 = quantile(obsfreqalert, 0.25, na.rm=TRUE),
    freqalert_Q2 = quantile(obsfreqalert, 0.5, na.rm=TRUE),
    freqalert_Q3 = quantile(obsfreqalert, 0.75, na.rm=TRUE),

    freqnotalert_Q1 = quantile(obsfreqnotalert, 0.25, na.rm=TRUE),
    freqnotalert_Q2 = quantile(obsfreqnotalert, 0.5, na.rm=TRUE),
    freqnotalert_Q3 = quantile(obsfreqnotalert, 0.75, na.rm=TRUE),

  )


write_rds(results1table, here::here("ord_2_2_3_results", "results1table.rds"))



bar_interval <-
  ggplot(data = results1table, aes(y = as.factor(tau), colour = as.factor(tau))) +

  geom_point(aes(x = interval_Q2), size = 2) +
  geom_errorbarh(aes(xmin = interval_Q1, xmax = interval_Q3), size = 1, height = 0.0) +
  #geom_text(aes(label = plyr::round_any(interval_Q1,0.01), x = interval_Q1), nudge_y = -.1) +
  #geom_text(aes(label = plyr::round_any(interval_Q2,0.01), x = interval_Q2), nudge_y = .1, vjust=0) +
  #geom_text(aes(label = plyr::round_any(interval_Q3,0.01), x = interval_Q3), nudge_y = -.1) +
  facet_grid(rows=vars(cluster)) +
  scale_y_discrete(breaks=c(0, 2^-(5:1)), labels = c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"))+
  scale_x_continuous(limits=c(0,1))+
  #scale_color_brewer(type='qual', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  labs(x = "Mean observation interval", y = "tau") +
  theme_bw()+
  theme(legend.position="none")

bar_freq <-
  ggplot(data = results1table, aes(y = as.factor(tau), colour = as.factor(tau))) +

  geom_point(aes(x = frequency_Q2), size = 2) +
  geom_errorbarh(aes(xmin = frequency_Q1, xmax = frequency_Q3), size = 1, height = 0.0) +
  #geom_text(aes(label = plyr::round_any(frequency_Q1,0.01), x = frequency_Q1), nudge_y = -.1) +
  #geom_text(aes(label = plyr::round_any(frequency_Q2,0.01), x = frequency_Q2), nudge_y = .1, vjust=0) +
  #geom_text(aes(label = plyr::round_any(frequency_Q3,0.01), x = frequency_Q3), nudge_y = -.1) +
  facet_grid(rows=vars(cluster)) +
  #scale_color_brewer(type='qual', direction=-1)+
  scale_colour_manual(values=viridisLite::viridis(n=6), na.value="transparent", labels=c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"), drop = FALSE)+
  scale_x_continuous(limits=c(0,4.1))+
  scale_y_discrete(breaks=c(0, 2^-(5:1)), labels = c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"))+
  labs(x = "Mean observation frequency (obs/day)", y = "tau") +
  theme_bw()+
  theme(legend.position="none")


bar_freqalert <-
  ggplot(data = results1table, aes(y = as.factor(tau), colour = as.factor(tau))) +

  geom_point(aes(x = freqalert_Q2), size = 2) +
  geom_errorbarh(aes(xmin = freqalert_Q1, xmax = freqalert_Q3), size = 1, height = 0.0) +
  #geom_text(aes(label = plyr::round_any(delay_Q1,0.01), x = delay_Q1), nudge_y = -.1) +
  #geom_text(aes(label = plyr::round_any(delay_Q2,0.01), x = delay_Q2), nudge_y = .1, vjust=0) +
  #geom_text(aes(label = plyr::round_any(delay_Q3,0.01), x = delay_Q3), nudge_y = -.1) +
  facet_grid(rows=vars(cluster)) +
  scale_color_brewer(type='qual', direction=-1)+
  scale_x_continuous(limits=c(0,4.1))+
  scale_y_discrete(breaks=c(0, 2^-(5:1)), labels = c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"))+
  labs(x = "Mean observation frequency in alert state (obs/day)", y = "tau") +
  theme_bw()+
  theme(legend.position="none")


bar_freqnotalert <-
  ggplot(data = results1table, aes(y = as.factor(tau), colour = as.factor(tau))) +

  geom_point(aes(x = freqnotalert_Q2), size = 2) +
  geom_errorbarh(aes(xmin = freqnotalert_Q1, xmax = freqnotalert_Q3), size = 1, height = 0.0) +
  #geom_text(aes(label = plyr::round_any(delay_Q1,0.01), x = delay_Q1), nudge_y = -.1) +
  #geom_text(aes(label = plyr::round_any(delay_Q2,0.01), x = delay_Q2), nudge_y = .1, vjust=0) +
  #geom_text(aes(label = plyr::round_any(delay_Q3,0.01), x = delay_Q3), nudge_y = -.1) +
  facet_grid(rows=vars(cluster)) +
  scale_color_brewer(type='qual', direction=-1)+
  scale_x_continuous(limits=c(0,4.1))+
  scale_y_discrete(breaks=c(0, 2^-(5:1)), labels = c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"))+
  labs(x = "Mean observation frequency not in alert state (obs/day)", y = "tau") +
  theme_bw()+
  theme(legend.position="none")


bar_delay <-
  ggplot(data = results1table, aes(y = as.factor(tau), colour = as.factor(tau))) +

  geom_point(aes(x = delay_Q2), size = 2) +
  geom_errorbarh(aes(xmin = delay_Q1, xmax = delay_Q3), size = 1, height = 0.0) +
  #geom_text(aes(label = plyr::round_any(delay_Q1,0.01), x = delay_Q1), nudge_y = -.1) +
  #geom_text(aes(label = plyr::round_any(delay_Q2,0.01), x = delay_Q2), nudge_y = .1, vjust=0) +
  #geom_text(aes(label = plyr::round_any(delay_Q3,0.01), x = delay_Q3), nudge_y = -.1) +
  facet_grid(rows=vars(cluster)) +
  scale_color_brewer(type='qual', direction=-1)+
  scale_x_continuous(limits=c(0,4.1))+
  scale_y_discrete(breaks=c(0, 2^-(5:1)), labels = c("0", "2^-5", "2^-4","2^-3","2^-2","2^-1"))+
  labs(x = "Mean delay before observing alert state (days)", y = "tau") +
  theme_bw()+
  theme(legend.position="none")

barplots <- cowplot::plot_grid(bar_freq, bar_freqalert, bar_freqnotalert, ncol=1, align="hv")



ggsave(here::here("ord_2_2_3_results", "barplots.png"), plot=barplots, width=20, height=20, units="cm")


# real data results -------------------------------------------------------


realresults <- read_rds(path=here::here("data", "simulations", "realresults_223.rds"))

realevaluate <- realresults %>%
  mutate(
    n=pmap_dbl(list(ptid, adp),
               function(ptid2, adp){
                  n <- nrow(adp)
                  #maxn <- nrow(filter(data_wide, ptid==ptid2))
                  ifelse(length(n)>0,n,NA_real_)
               }
               ),
    maxn=pmap_dbl(list(ptid, adp),
                  function(ptid2, adp){
                    dat <- filter(data_wide,
                                  ptid == ptid2,
                                  fup_days > min(adp$time_J)
                                  )
                    nrow(dat)
                  }
                  ),
    prop=n/maxn
  ) %>%
  select(-adp)



realevaluate %>%
  ggplot()+
  facet_grid(rows=(vars(tau)))+
  geom_density(aes(x=prop), alpha=0.6)+
  theme_bw()

realevaluate %>%
  ggplot()+
  geom_jitter(aes(x=prop, y=as.factor(tau)), width=0, height=0.2, alpha=0.6)+
  theme_bw()


plot_used <- realevaluate %>%
  ggplot()+
  geom_line(aes(x=tau, y=prop+runif(length(prop),-0.02,0), group=ptid), alpha=0.2)+
  labs(x="tau", y="Proportion of observations used")+
  theme_bw()


ggsave(here::here("ord_2_2_3_results", "plot_obsused.png"), plot=plot_used, width=20, height=20, units="cm")


