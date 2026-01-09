library(tidyverse)
library(mipfp) # package for raking
library(rstanarm)

# As weighted mean of the values in Table 1 is all for simple poststratification, 
# no codes are provided for simple poststratification. 


# Raking ------------------------------------------------------------------

# We first preprocess KNHANES marginal count data, previously obtained from SAS codes.

# 1. KNHANES bootstrapped marginal data preprocessing --------------------------------------------------------

# As marginal counts are not integer, their sums may not strictly equal.
# e.g. sum of marginal counts by age groups may differ from that of sex groups.
# Before applying raking, we should make sums equal. 

reg2030 = read_csv("/knhanes/knhanes_boot_marginal_2030.csv")
reg40 = read_csv("/knhanes/knhanes_boot_marginal_40.csv")

# 1. For each boot_iteration, check whether rounded sums of marginal counts by variables are equal.

dup_2030 = reg2030 %>% 
  group_by(boot_iteration, label) %>% 
  summarise(n = round(sum(count))) %>% 
  distinct(n) %>% pull(boot_iteration)

dup_40 = reg40 %>% 
  group_by(boot_iteration, label) %>% 
  summarise(n = round(sum(count))) %>% 
  distinct(n) %>% pull(boot_iteration)

# Among 1000 bootstrapped samples, sex and smk of sample 301 had different rounded value. 
# modified by adding 0.1
reg40[reg40$boot_iteration==301 & reg40$label=='sex', 'count'] = reg40[reg40$boot_iteration==301 & reg40$label=='sex', 'count'] + 0.1
reg40[reg40$boot_iteration==301 & reg40$label=='smk', 'count'] = reg40[reg40$boot_iteration==301 & reg40$label=='smk', 'count'] + 0.1


# 2. Modify every counts into integer, to make sums equal without rounding. 

largest_remainder_round = function(x) {
  # Calculate the target sum (rounded)
  target_sum = round(sum(x))
  
  # Get integer parts
  int_parts = floor(x)
  
  # Get fractional parts
  frac_parts = x - int_parts
  
  # Calculate how many values need to be rounded up
  n_round_up = target_sum - sum(int_parts)
  
  # Round up the values with largest fractional parts
  result = int_parts
  if (n_round_up > 0) {
    indices_to_round_up = order(frac_parts, decreasing = TRUE)[1:n_round_up]
    result[indices_to_round_up] = result[indices_to_round_up] + 1
  }
  
  return(as.integer(result))
}

reg2030_rounded = reg2030 %>%
  group_by(boot_iteration, label) %>%
  mutate(count_int = largest_remainder_round(count)) %>%
  ungroup()

reg40_rounded = reg40 %>%
  group_by(boot_iteration, label) %>%
  mutate(count_int = largest_remainder_round(count)) %>%
  ungroup()


# 3. As some region counts are not presented in data due to 0 count, make that row. 
# check the region which do not have 1000 rows.

table(reg2030[reg2030$label=='reg','cat'])
table(reg40[reg40$label=='reg','cat'])

# -> only region 7 has some missings. 

missing7_2030 = setdiff(1:1000, reg2030 %>% filter(label=='reg' & cat== 7) %>% pull(boot_iteration))
missing7_40 = setdiff(1:1000, reg40 %>% filter(label=='reg' & cat== 7) %>% pull(boot_iteration))
temp_2030 = data.frame(boot_iteration = missing7_2030, label = 'reg', cat = 7, count = 0, count_int = 0)
temp_40 = data.frame(boot_iteration = missing7_40, label = 'reg', cat = 7, count = 0, count_int = 0)

reg2030_final = rbind(reg2030_rounded, temp_2030)
reg40_final = rbind(reg40_rounded, temp_40)

reg2030_final = reg2030_final %>% 
  mutate(label = factor(label, levels = c("age", 'sex', 'smk', 'drk', 'reg'))) %>%
  arrange(boot_iteration, label, cat) %>%
  mutate(label = as.character(label))

reg40_final = reg40_final %>% 
  mutate(label = factor(label, levels = c("age", 'sex', 'smk', 'drk', 'reg'))) %>%
  arrange(boot_iteration, label, cat) %>%
  mutate(label = as.character(label))


write.csv(reg2030_final, "/knhanes/knhanes_boot_marginal_2030_augmented.csv")
write.csv(reg40_final, "/knhanes/knhanes_boot_marginal_40_augmented.csv")

# do the same preprocessing for original KNHANES(not bootstrapped) marginal counts

# 2. raking ------------------------------------------------------------------

# Now we conduct raking using preprocessed data. 

# assume we have extracted cellwise count and obesity data of nhis-nsc 1000 bootstrap samples.

nhis_2030 = read.csv("/nhis_nsc/nhis_boot_cellwise_obesity_2030.csv")
nhis_40 = read.csv("/nhis_nsc/nhis_boot_cellwise_obesity_40.csv")

# 1. Get a raking point estimate 

# import marginal counts of original KNHANES
marginal_knhanes_2030 = read.csv("/knhanes/knhanes_marginal_2030.csv")
marginal_knhanes_40 = read.csv("/knhanes/knhanes_marginal_40.csv")

# 1-1 conduct raking with original nhis-nsc data

dims_40 = c(5,2,2,2) # in the order of age, sex, smoke, drink
dims_2030 = c(4,2,2,2) 

target_data_2030 = list(marginal_knhanes_2030 %>% filter(label=='age') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='sex') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='smk') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='drk') %>% pull(count)
)

target_data_40 = list(marginal_knhanes_40 %>% filter(label=='age') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='sex') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='smk') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='drk') %>% pull(count)
)

# boot_iteration == 0 represents the cellwise counts of original nhis-nsc data
sample_count_wo_region_2030 = nhis_2030 %>% filter(boot_iteration ==0) %>% arrange(drink, smoke, SEX, age_cat) %>% pull(count)
sample_count_wo_region_40 = nhis_40 %>% filter(boot_iteration ==0) %>% arrange(drink, smoke, SEX, age_cat) %>% pull(count)

# raking starts from nhis-nsc counts
seed_array_2030 = array(sample_count_wo_region_2030, dim=dims_2030)
seed_array_40 = array(sample_count_wo_region_40, dim=dims_40)


ipf_result_2030 = Ipfp(seed = seed_array_2030,
                  target.list = list(1, 2, 3, 4),
                  target.data = target_data_2030,
                  print = FALSE)

ipf_result_40 = Ipfp(seed = seed_array_40,
                       target.list = list(1, 2, 3, 4),
                       target.data = target_data_40,
                       print = FALSE)

# cellwise population counts estimated through raking 
raking_2030 = array(ipf_result_2030$x.hat, dim = 32)
raking_40 = array(ipf_result_40$x.hat, dim = 40)

# 1-2 apply poststratification 
nhis_obesity_2030 = nhis_2030 %>% filter(boot_iteration ==0) %>% 
  arrange(drink, smoke, SEX, age_cat) %>% pull(obesity)

nhis_obesity_40 = nhis_40 %>% filter(boot_iteration ==0) %>% 
  arrange(drink, smoke, SEX, age_cat) %>% pull(obesity)

weighted.mean(nhis_obesity_2030, raking_2030)
weighted.mean(nhis_obesity_40, raking_40)


# 2. Get bootstrap CIs

# import previously preprocessed knhanes bootstrap marginal count data

knhanes_boot_marginal_2030 = read.csv("/knhanes/knhanes_boot_marginal_2030_augmented.csv")
knhanes_boot_marginal_40 = read.csv("/knhanes/knhanes_boot_marginal_40_augmented.csv")

# Repeat above process for each bootstrap samples
raking_boot_2030 = c()
raking_boot_40 = c()

for (i in 1:1000){
  temp_2030 = knhanes_boot_marginal_2030[knhanes_boot_marginal_2030$boot_iteration==i,]
  temp_40 = knhanes_boot_marginal_40[knhanes_boot_marginal_40$boot_iteration==i,]
  
  target_data_2030 = list(temp_2030 %>% filter(label=='age') %>% pull(count_int),
                          temp_2030 %>% filter(label=='sex') %>% pull(count_int),
                          temp_2030 %>% filter(label=='smk') %>% pull(count_int),
                          temp_2030 %>% filter(label=='drk') %>% pull(count_int)
  )
  
  target_data_40 = list(temp_40 %>% filter(label=='age') %>% pull(count_int),
                        temp_40 %>% filter(label=='sex') %>% pull(count_int),
                        temp_40 %>% filter(label=='smk') %>% pull(count_int),
                        temp_40 %>% filter(label=='drk') %>% pull(count_int)
  )
  
  nhis_count_2030 = nhis_2030 %>% filter(boot_iteration ==i) %>% 
    arrange(drink, smoke, SEX, age_cat) %>% pull(count)
  nhis_count_40 = nhis_40 %>% filter(boot_iteration ==i) %>% 
    arrange(drink, smoke, SEX, age_cat) %>% pull(count)
  
  seed_array_2030 = array(nhis_count_2030, dim=dims_2030)
  seed_array_40 = array(nhis_count_40, dim=dims_40)
  
  ipf_result_2030 = Ipfp(seed = seed_array_2030,
                         target.list = list(1, 2, 3, 4),
                         target.data = target_data_2030,
                         print = FALSE)
  
  ipf_result_40 = Ipfp(seed = seed_array_40,
                       target.list = list(1, 2, 3, 4),
                       target.data = target_data_40,
                       print = FALSE)
  
  raking_2030 = array(ipf_result_2030$x.hat, dim = 32)
  raking_40 = array(ipf_result_40$x.hat, dim = 40)
  
  nhis_obesity_2030 = nhis_2030 %>% filter(boot_iteration ==i) %>% 
    arrange(drink, smoke, SEX, age_cat) %>% pull(obesity)
  nhis_obesity_40 = nhis_40 %>% filter(boot_iteration ==i) %>% 
    arrange(drink, smoke, SEX, age_cat) %>% pull(obesity)
  
  raking_boot_2030 = c(raking_boot_2030, weighted.mean(nhis_obesity_2030, raking_2030))
  raking_boot_40 = c(raking_boot_40, weighted.mean(nhis_obesity_40, raking_40))
}

raking_results = data.frame(cbind(raking_boot_2030, raking_boot_40))

quantile(raking_results$raking_boot_2030, c(0.025, 0.975))
quantile(raking_results$raking_boot_40, c(0.025, 0.975))


# MRP ---------------------------------------------------------------------


# 1. NHIS-NSC data preprocessing ---------------------------------------------

# import nhis-nsc data
nhis_2030 = read.csv("/nhis_nsc/nhis_boot_cellwise_obesity_2030.csv")
nhis_40 = read.csv("/nhis_nsc/nhis_boot_cellwise_obesity_40.csv")

# As the region code differed between KNHANES and NHIS-NSC, modify them. 
nhis_2030 = nhis_2030 %>% 
  mutate(region2 = case_when(
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
  )) %>% 
  filter(boot_iteration == 0) # choose only original nhis-nsc data

nhis_40 = nhis_40 %>% 
  mutate(region2 = case_when(
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
  )) %>% 
  filter(boot_iteration == 0) 


# 2. Bayesian multilevel regression ----------------------------------------------------------------

# fit bayesian multilevel regression model
set.seed(1)
fit_2030 = stan_glmer(obese ~ factor(age_cat) + factor(SEX) + factor(smoke) + factor(drink) + 
                        (1|region2), 
                      family = binomial, 
                      weights = count, 
                      data = nhis_2030, 
                      cores = 4)

set.seed(1)
fit_40 = stan_glmer(obese ~ factor(age_cat) + factor(SEX) + factor(smoke) + factor(drink) + 
                      (1|region2), 
                    family = binomial, 
                    weights = count, 
                    data = nhis_40, 
                    cores = 4)

# predict cellwise obesity prevalence
pred_matrix_2030 = posterior_epred(
  fit_2030, 
  newdata = expand.grid(age_cat = c(2,2.5,3,3.5), SEX = 1:2, smoke = 0:1, 
                        drink = 0:1, region2 = 1:17),
  re.form = NULL  
)

pred_matrix_40 = posterior_epred(
  fit_40, 
  newdata = expand.grid(age_cat = c(4,5,6,7,8), SEX = 1:2, smoke = 0:1, 
                        drink = 0:1, region2 = 1:17),
  re.form = NULL  
)


# keep only last 250 estimates from each chain.
pred_2030 = data.frame(pred_matrix_2030)[c(751:1000, 1751:2000, 2751:3000, 3751:4000),]
pred_40 = data.frame(pred_matrix_40)[c(751:1000, 1751:2000, 2751:3000, 3751:4000),]


# 3. Poststratification - point estimate ----------------------------------------------------------

# import original KNHANES marginal counts
marginal_knhanes_2030 = read.csv("/knhanes/knhanes_marginal_2030.csv")
marginal_knhanes_40 = read.csv("/knhanes/knhanes_marginal_40.csv")

target_data_2030 = list(marginal_knhanes_2030 %>% filter(label=='age') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='sex') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='smk') %>% pull(count),
                        marginal_knhanes_2030 %>% filter(label=='drk') %>% pull(count), 
                        marginal_knhanes_2030 %>% filter(label=='reg') %>% pull(count) 
)

target_data_40 = list(marginal_knhanes_40 %>% filter(label=='age') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='sex') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='smk') %>% pull(count),
                      marginal_knhanes_40 %>% filter(label=='drk') %>% pull(count), 
                      marginal_knhanes_40 %>% filter(label=='reg') %>% pull(count) 
)

# raking start from uniform distribution
seed_array_2030 = array(1, dim=dims_2030)
seed_array_40 = array(1, dim=dims_40)

ipf_result_2030 = Ipfp(seed = seed_array_2030,
                       target.list = list(1, 2, 3, 4, 5),
                       target.data = target_data_2030,
                       print = FALSE)

ipf_result_40 = Ipfp(seed = seed_array_40,
                     target.list = list(1, 2, 3, 4, 5),
                     target.data = target_data_40,
                     print = FALSE)

raking_2030 = array(ipf_result_2030$x.hat, dim = 4*2*2*2*17)
raking_40 = array(ipf_result_40$x.hat, dim = 5*2*2*2*17)

point_2030 = c()
point_40 = c()

# poststratify previously obtained cellwise obesity estimates 1000 times
for (i in 1:1000){
  b = pred_2030[i, ]
  point_2030 = c(point_2030, weighted.mean(b, raking_2030) )
  
  b = pred_40[i, ]
  point_40 = c(point_40, weighted.mean(b, raking_40) )
}

# get the MRP point estimate as a mean
mean(point_2030) 
mean(point_40)



# 4. Poststratification - bootstrap CIs ----------------------------------------------------------

# import KNHANES bootstrapped marginal counts
knhanes_boot_marginal_2030 = read.csv("/knhanes/knhanes_boot_marginal_2030_augmented.csv")
knhanes_boot_marginal_40 = read.csv("/knhanes/knhanes_boot_marginal_40_augmented.csv")

raking_boot_2030 = rep(NA, 4*2*2*2*17)
raking_boot_40 = rep(NA, 5*2*2*2*17)

dims_2030 = c(4,2,2,2,17)
dims_40 = c(5,2,2,2,17)

# iterate raking 1000 times 
for (i in 1:1000){
  temp_2030 = knhanes_boot_marginal_2030[knhanes_boot_marginal_2030$boot_iteration==i,]
  temp_40 = knhanes_boot_marginal_40[knhanes_boot_marginal_40$boot_iteration==i,]
  
  target_data_2030 = list(temp_2030 %>% filter(label=='age') %>% pull(count_int),
                          temp_2030 %>% filter(label=='sex') %>% pull(count_int),
                          temp_2030 %>% filter(label=='smk') %>% pull(count_int),
                          temp_2030 %>% filter(label=='drk') %>% pull(count_int), 
                          temp_2030 %>% filter(label=='reg') %>% pull(count_int) 
  )
  
  target_data_40 = list(temp_40 %>% filter(label=='age') %>% pull(count_int),
                        temp_40 %>% filter(label=='sex') %>% pull(count_int),
                        temp_40 %>% filter(label=='smk') %>% pull(count_int),
                        temp_40 %>% filter(label=='drk') %>% pull(count_int), 
                        temp_40 %>% filter(label=='reg') %>% pull(count_int) 
  )
  
  # start from uniform distribution
  seed_array_2030 = array(1, dim=dims_2030)
  seed_array_40 = array(1, dim=dims_40)
  
  ipf_result_2030 = Ipfp(seed = seed_array_2030,
                         target.list = list(1, 2, 3, 4, 5),
                         target.data = target_data_2030,
                         print = FALSE)
  
  ipf_result_40 = Ipfp(seed = seed_array_40,
                       target.list = list(1, 2, 3, 4, 5),
                       target.data = target_data_40,
                       print = FALSE)
  
  raking_2030 = array(ipf_result_2030$x.hat, dim = 4*2*2*2*17)
  raking_40 = array(ipf_result_40$x.hat, dim = 5*2*2*2*17)
  
  raking_boot_2030 = rbind(raking_boot_2030, raking_2030)
  raking_boot_40 = rbind(raking_boot_40, raking_40)
}

raking_boot_2030 = data.frame(raking_boot_2030[2:1001,])
raking_boot_40 = data.frame(raking_boot_40[2:1001,])

mean_2030 = c()
mean_40 = c()

# Get CIs of point estimate 
for (i in 1:1000){
  a = raking_boot_2030[i, ]
  b = pred_2030[i, ]
  mean_2030 = c(mean_2030, weighted.mean(b, a) )
  
  a = raking_boot_40[i, ]
  b = pred_40[i, ]
  mean_40 = c(mean_40, weighted.mean(b, a) )
  
}

quantile(mean_2030, c(0.025, 0.975)) 
quantile(mean_40, c(0.025, 0.975)) 
