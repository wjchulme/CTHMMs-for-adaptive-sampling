###############

# all previous data-processing steps are hidden here to protect participant privacy (includes participation dates, medication, sex, etc..)
#all data needed for this project is contained in the data_canon data.frame
#write_rds(data_canon, path = here::here("data","data_canon.rds"))

###############

read_rds(path = here::here("data","data_canon.rds"))

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


write_rds(data_wide, path = here::here("data","data_wide.rds"))








# explore -----------------------------------------------------------------




# how many times was each question asked+answered per patient over the course of the study?

expl_question_count <- data_canon %>%
  count(ptid, itemgrp)

# when were questions answered?
data_canon %>%
  ggplot() +
  geom_jitter(aes(x=fup_days%%1, y=factor(itemgrp)), height=0.1, size=0.5, alpha=0.5)+
  scale_x_time()+
  facet_wrap(ptid~.)+
  theme_bw()+
  theme(legend.position="bottom", legend.direction="horizontal", strip.text = element_blank())+
  NULL

data_canon %>%
  ggplot() +
  geom_point(aes(x=fup_days, y=fup_days%%1, colour=score), size=0.1)+
  facet_wrap(~ptid)+
  facet_grid(itemgrp~.)+
  scale_y_time()+
  theme_bw()+
  NULL

data_canon %>%
  group_by(ptid, itemgrp) %>%
  mutate(timediff=fup_days-lag(fup_days)) %>%
  ggplot() +
  geom_histogram(aes(x=timediff/(60*60)))+
  NULL

# responses per day -------------------------------------------------------



# absolute scores per patient over time

data_canon %>%
  ggplot() +
  geom_jitter(aes(x=fup_days, y=score, colour=itemgrp), width=0, height=0.1, alpha=0.1, size=0.8) +
  #geom_point(aes(x=fup_days, y=score, colour=itemgrp), alpha=0.2, size=0.8) +
  facet_wrap(~ptid)+
  scale_colour_brewer(palette="Dark2")+ # or Set2
  ylim(c(0.8,7.2))+
  theme_bw()+
  theme(legend.position="bottom",legend.direction="horizontal", strip.text = element_blank())



fig_score_encounter <-
data_canon %>%
ggplot() +
  geom_tile(aes(x=fup_days, y=itemgrp, fill=as.factor(score)), height=0.6, width=0.5, colour="transparent") +
  facet_wrap(~ptid, ncol=10)+
  scale_fill_manual(values=c(viridisLite::viridis(n=7)), na.value="transparent", breaks=1:7)+
  labs(x="Day", y=NULL, fill="Response")+
  theme_bw()+
  theme(
    legend.position="bottom",legend.direction="horizontal",
    strip.text = element_blank(),
    panel.grid.major.y = element_blank()
  )+
  guides(fill=guide_legend(nrow=1))
ggsave(filename="figures/raw trajectories.png", fig_score_encounter, units='cm', height=20, width=30)

