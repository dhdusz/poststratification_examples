library(tidyverse)
library(survey)
library(readxl)
library(haven)

# In implementing IPSW in obesity example, we need counts of variable combination * obesity status. 
# i.e. sample / population count of each age * sex * smk * drk * region * obesity combination.
# suppose we have that nhis-nsc data

nhis_2030 = read.csv("/nhis_nsc/nhis_cellwise_obese_counts_2030.csv")
nhis_40 = read.csv("/nhis_nsc/nhis_cellwise_obese_counts_40.csv")  
  

# Modify nhis-nsc data to match complex survey structure of KNHANES.

# set domain2030 or domain40 accordingly (cf. SAS codes)
# designate different psu's for every different individual
# designate single stratum for every participants
# weight = 1 for all individual
# Ri = 1 ("participated") for all nhis-nsc participants 

nhis_2030 = nhis_2030 %>% uncount(count) %>% mutate(domain2030 = 1,
                                                    domain40 = 0,
                                                    psu = paste0('Z',row_number()), 
                                                    kstrata = 1000, 
                                                    wt_itvex = 1, 
                                                    Ri = 1) %>% rename(sex = SEX)
nhis_40 = nhis_40 %>% uncount(count) %>% mutate(domain2030 = 0,
                                                domain40 = 1,
                                                psu = paste0('Y',row_number()), 
                                                kstrata = 1001, 
                                                wt_itvex = 1, 
                                                Ri = 1) %>% rename(sex = SEX)

nhis_2030 = nhis_2030 %>% 
  mutate(region = case_when(
    region == 11 ~ 1,
    region == 26 ~ 2,
    region == 27 ~ 3,
    region == 28 ~ 4,
    region == 29 ~ 5,
    region == 30 ~ 6,
    region == 31 ~ 7,
    region == 36 ~ 8,
    region == 41 ~ 9,
    region == 42 ~ 10,
    region == 43 ~ 11,
    region == 44 ~ 12,
    region == 45 ~ 13,
    region == 46 ~ 14,
    region == 47 ~ 15,
    region == 48 ~ 16,
    region == 50 ~ 17
  ))

nhis_40 = nhis_40 %>% 
  mutate(region = case_when(
    region == 11 ~ 1,
    region == 26 ~ 2,
    region == 27 ~ 3,
    region == 28 ~ 4,
    region == 29 ~ 5,
    region == 30 ~ 6,
    region == 31 ~ 7,
    region == 36 ~ 8,
    region == 41 ~ 9,
    region == 42 ~ 10,
    region == 43 ~ 11,
    region == 44 ~ 12,
    region == 45 ~ 13,
    region == 46 ~ 14,
    region == 47 ~ 15,
    region == 48 ~ 16,
    region == 50 ~ 17
  ))

# For KNHANES, set Ri = 0, and delete those with weight = NA

knhanes = read.csv("/knhanes/knhanes19.csv")

knhanes = knhanes %>% select(-ID) %>% mutate(Ri = 0) %>% filter(!is.na(wt_itvex))

# generate aggregate data
aggregate2030 = rbind(knhanes, nhis_2030)
aggregate40 = rbind(knhanes, nhis_40)

# generate svydesign
design_full_2030 = 
  svydesign(
  id = ~psu, 
  strata = ~kstrata, 
  weights = ~wt_itvex, 
  data = aggregate2030, 
  nest = TRUE
)

design_full_40 = 
  svydesign(
    id = ~psu, 
    strata = ~kstrata, 
    weights = ~wt_itvex, 
    data = aggregate40, 
    nest = TRUE
  )

# generate domain 
design_2030 = subset(design_full_2030, domain2030 ==1)
design_40 = subset(design_full_40, domain40 ==1)


# fit svyglm using different sets of variables 
svytry2030 = svyglm(Ri~factor(age_cat) + sex, 
                    design = design_2030, family = quasibinomial())
svytry2030 = svyglm(Ri~factor(age_cat) + sex + smoke + drink, 
                    design = design_2030, family = quasibinomial())
svytry2030 = svyglm(Ri~factor(age_cat) + sex + smoke + drink + factor(region), 
                    design = design_2030, family = quasibinomial())

svytry40 = svyglm(Ri~factor(age_cat) + sex, 
                  design = design_40, family = quasibinomial())
svytry40 = svyglm(Ri~factor(age_cat) + sex + smoke + drink, 
                  design = design_40, family = quasibinomial())
svytry40 = svyglm(Ri~factor(age_cat) + sex + smoke + drink + factor(region), 
                  design = design_40, family = quasibinomial())

# prediction
predict2030_nhis = as.numeric(predict(svytry2030, newdata = nhis_2030, type = 'response'))
predict40_nhis = as.numeric(predict(svytry40, newdata = nhis_40, type = 'response'))

predict2030_knhanes = as.numeric(predict(svytry2030, 
                                         newdata = knhanes %>% filter(domain2030==1), 
                                         type = 'response'))
predict40_knhanes = as.numeric(predict(svytry40, 
                                       newdata = knhanes %>% filter(domain40==1), 
                                       type = 'response'))

# Get mu_ALP estimate 
nhis_2030$pi = predict2030_nhis
nhis_40$pi = predict40_nhis

nhis_2030$wi = (1-nhis_2030$pi) / nhis_2030$pi
nhis_40$wi =(1- nhis_40$pi) / nhis_40$pi

knhanes[knhanes$domain2030==1,'pi_2030'] = predict2030_knhanes
knhanes[knhanes$domain40==1,'pi_40'] = predict40_knhanes

knhanes$wi_2030 = (1-knhanes$pi_2030) /  knhanes$pi_2030
knhanes$wi_40 = (1-knhanes$pi_40) /  knhanes$pi_40

mu_ALP_2030 = weighted.mean(nhis_2030$obese, nhis_2030$wi)
mu_ALP_40 = weighted.mean(nhis_40$obese, nhis_40$wi)


# variance estimation - bootstrap -----------------------------------------

# As variance estimation suggested by Wang et al. resulted in too small variance, 
# we instead applied bootstrap variance. 

results_mu = data.frame(matrix(ncol = 6, nrow = 1000))

set.seed(1)
for (i in 1:1000){
  print(i)
  # generate bootstrap sample
  boot_ind_2030 = sample(1:nrow(nhis_2030), replace=TRUE)
  boot_ind_40 = sample(1:nrow(nhis_40), replace=TRUE)
  
  nhis_boot_2030 = nhis_2030[boot_ind_2030,]
  nhis_boot_40 = nhis_40[boot_ind_40,]
  
  knhanes_boot = read_sas(paste0("/knhanes_boot/boot_",i,".sas7bdat"))
  knhanes_boot = knhanes_boot %>% 
    filter(!is.na(wt_itvex)) %>% 
    mutate(wt_itvex = bootstrap_weight, 
           Ri = 0) %>%
    select(age_cat, sex, smoke, drink, region, obese, domain2030, domain40, psu, kstrata, wt_itvex, Ri)
                          
  aggregate2030 = rbind(knhanes_boot, nhis_boot_2030)
  aggregate40 = rbind(knhanes_boot, nhis_boot_40)
  
  # repeat above procedures

  design_full_2030 = 
    svydesign(
      id = ~psu, 
      strata = ~kstrata, 
      weights = ~wt_itvex, 
      data = aggregate2030, 
      nest = TRUE
    )
  
  design_full_40 = 
    svydesign(
      id = ~psu, 
      strata = ~kstrata, 
      weights = ~wt_itvex, 
      data = aggregate40, 
      nest = TRUE
    )
  
  # domain 생성 
  design_2030 = subset(design_full_2030, domain2030 ==1)
  design_40 = subset(design_full_40, domain40 ==1)
  
  # svyglm fitting 
  svytry2030_as = svyglm(Ri~factor(age_cat) + sex, 
                      design = design_2030, family = quasibinomial())
  svytry2030_assd = svyglm(Ri~factor(age_cat) + sex + smoke + drink, 
                      design = design_2030, family = quasibinomial())
  svytry2030_assdr = svyglm(Ri~factor(age_cat) + sex + smoke + drink + factor(region), 
                      design = design_2030, family = quasibinomial())
  
  svytry40_as = svyglm(Ri~factor(age_cat) + sex, 
                    design = design_40, family = quasibinomial())
  svytry40_assd = svyglm(Ri~factor(age_cat) + sex + smoke + drink, 
                    design = design_40, family = quasibinomial())
  svytry40_assdr = svyglm(Ri~factor(age_cat) + sex + smoke + drink + factor(region), 
                    design = design_40, family = quasibinomial())
  
  # prediction
  predict2030_nhis_as = as.numeric(predict(svytry2030_as, newdata = nhis_boot_2030, type = 'response'))
  predict2030_nhis_assd = as.numeric(predict(svytry2030_assd, newdata = nhis_boot_2030, type = 'response'))
  predict2030_nhis_assdr = as.numeric(predict(svytry2030_assdr, newdata = nhis_boot_2030, type = 'response'))
  
  predict40_nhis_as = as.numeric(predict(svytry40_as, newdata = nhis_boot_40, type = 'response'))
  predict40_nhis_assd = as.numeric(predict(svytry40_assd, newdata = nhis_boot_40, type = 'response'))
  predict40_nhis_assdr = as.numeric(predict(svytry40_assdr, newdata = nhis_boot_40, type = 'response'))
  
  predict2030_knhanes_as = as.numeric(predict(svytry2030_as, 
                                           newdata = knhanes_boot %>% filter(domain2030==1), 
                                           type = 'response'))
  predict2030_knhanes_assd = as.numeric(predict(svytry2030_assd, 
                                              newdata = knhanes_boot %>% filter(domain2030==1), 
                                              type = 'response'))
  predict2030_knhanes_assdr = as.numeric(predict(svytry2030_assdr, 
                                              newdata = knhanes_boot %>% filter(domain2030==1), 
                                              type = 'response'))
  
  predict40_knhanes_as = as.numeric(predict(svytry40_as, 
                                         newdata = knhanes_boot %>% filter(domain40==1), 
                                         type = 'response'))
  predict40_knhanes_assd = as.numeric(predict(svytry40_assd, 
                                            newdata = knhanes_boot %>% filter(domain40==1), 
                                            type = 'response'))
  predict40_knhanes_assdr = as.numeric(predict(svytry40_assdr, 
                                            newdata = knhanes_boot %>% filter(domain40==1), 
                                            type = 'response'))
  
  nhis_boot_2030$pi_as = predict2030_nhis_as
  nhis_boot_2030$pi_assd = predict2030_nhis_assd
  nhis_boot_2030$pi_assdr = predict2030_nhis_assdr
  
  nhis_boot_40$pi_as = predict40_nhis_as
  nhis_boot_40$pi_assd = predict40_nhis_assd
  nhis_boot_40$pi_assdr = predict40_nhis_assdr
  
  nhis_boot_2030$wi_as = (1-nhis_boot_2030$pi_as) / nhis_boot_2030$pi_as
  nhis_boot_2030$wi_assd = (1-nhis_boot_2030$pi_assd) / nhis_boot_2030$pi_assd
  nhis_boot_2030$wi_assdr = (1-nhis_boot_2030$pi_assdr) / nhis_boot_2030$pi_assdr
  
  nhis_boot_40$wi_as =(1- nhis_boot_40$pi_as) / nhis_boot_40$pi_as
  nhis_boot_40$wi_assd =(1- nhis_boot_40$pi_assd) / nhis_boot_40$pi_assd
  nhis_boot_40$wi_assdr =(1- nhis_boot_40$pi_assdr) / nhis_boot_40$pi_assdr
  
  mu_ALP_2030_as = weighted.mean(nhis_boot_2030$obese, nhis_boot_2030$wi_as)
  mu_ALP_2030_assd = weighted.mean(nhis_boot_2030$obese, nhis_boot_2030$wi_assd)
  mu_ALP_2030_assdr = weighted.mean(nhis_boot_2030$obese, nhis_boot_2030$wi_assdr)
  
  mu_ALP_40_as = weighted.mean(nhis_boot_40$obese, nhis_boot_40$wi_as)
  mu_ALP_40_assd = weighted.mean(nhis_boot_40$obese, nhis_boot_40$wi_assd)
  mu_ALP_40_assdr = weighted.mean(nhis_boot_40$obese, nhis_boot_40$wi_assdr)

  mu_s = c(mu_ALP_2030_as, mu_ALP_2030_assd, mu_ALP_2030_assdr, mu_ALP_40_as, mu_ALP_40_assd, mu_ALP_40_assdr)
  
  results_mu[i, ] = mu_s
  
}

write.csv(results_mu, "C:\\Users\\fermat\\OneDrive\\바탕 화면\\과제\\과제\\poststratification\\ipsw_boot_results.csv")

quantile(results_mu$X1, c(0.025, 0.975))
quantile(results_mu$X2, c(0.025, 0.975))
quantile(results_mu$X3, c(0.025, 0.975))
quantile(results_mu$X4, c(0.025, 0.975))
quantile(results_mu$X5, c(0.025, 0.975))
quantile(results_mu$X6, c(0.025, 0.975))
