
library(Hmisc)
library(tidyverse)
library(haven)


# Load deaths from VR countries and prepare dataset

causes <- c("preterm", "intrapartum", "congenital", "sepsis", "pneumonia", "injuries", "other")

death0 <- read_dta("DATA/good_VRneonates_observed_20250828.dta")

deaths <- transmute(death0, iso3=isocode, year, country, period_vr, preterm, intrapartum=asphyxia, congenital, sepsis=othinfection, pneumonia=ari, injuries, other=(other1+other2+hiv+malaria+pertussis+meningitis+tetanus)) %>%
  arrange(iso3, year, period_vr) %>% mutate(recnr=1:n()) %>%
  pivot_longer(cols=preterm:other, names_to="cause", values_to="n") %>%
  mutate(cause=factor(cause, levels=causes),
         preterm=ifelse(cause=="preterm",1,0),
         intrapartum=ifelse(cause=="intrapartum",1,0),
         congenital=ifelse(cause=="congenital",1,0),
         sepsis=ifelse(cause=="sepsis",1,0),
         pneumonia=ifelse(cause=="pneumonia",1,0),
         injuries=ifelse(cause=="injuries",1,0),
         other=ifelse(cause=="other",1,0)) %>%
  arrange(period_vr, iso3, year, cause)

deaths_e <- filter(deaths, period_vr==2) %>% mutate(row=c(1:n()))

deaths_l <- filter(deaths, period_vr==4) %>% mutate(row=c(1:n()))



# Load covariates data from VR countries and prepare dataset

covar0 <- read_csv("DATA/CovariateDatabase2023-wide_20250728.csv")

studies <- transmute(covar0, iso3, year, country=country_name, 
                     intercept = 1,
                     gni=gni_sm,
                     gfr=ifelse(gfr_sm<0.001, gfr_sm*1000, gfr_sm),
                     gini=gini_sm,
                     u5mr=u5mr_sm, nmr=nmr_sm, lbw=lbw_sm*100,
                     dpt=vac_dtp3_sm*100, 
                     anc=anc4_sm*100, 
                     femlit=literacy_f_sm*100) 

studies_e <- group_by(deaths_e, iso3, year, recnr) %>%
  summarise(first=min(row), last=max(row), N=sum(n)) %>% right_join(studies,.) %>%
  mutate(cv_group=rep(sample(c(1:10)), length.out=nrow(.)))

studies_l <- group_by(deaths_l, iso3, year, recnr) %>%
  summarise(first=min(row), last=max(row), N=sum(n)) %>% right_join(studies,.) %>%
  mutate(cv_group=rep(sample(c(1:10)), length.out=nrow(.)))



# Covariables names and labels

vdt <- c("preterm", "intrapartum", "congenital", "sepsis", "pneumonia", "injuries", "other")
vxn <- c("gni", "gfr", "gini", "u5mr", "nmr", "lbw", "dpt", "anc", "femlit")
vxr <- c("iso3") 
lxf <- c("GNI per capita", "General fertility rate", "Gini index","Under-five mortality rate", "Neonatal mortality rate", "Low birthweight rate", "Diphtheria/pertussis/tetanus vaccine coverage", "Antenatal care coverage (at least 1 visit)", "Adult female literacy rate")


save(studies_e, studies_l, deaths_e, deaths_l, vdt, vxn, vxr, lxf, file="DATA/data_neonates_LM_20250830.RData")


