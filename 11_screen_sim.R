library(pgxsim)
library(tidyverse)
my_seed <- 10001
set.seed(my_seed)

#define the parameters that will remain fixed
fixed_df <- data_frame(type='discrete', mu=.5, lb=-4, ub=Inf, minconc=0.001, maxconc=30,
                       ndoses=c(7,10), nreps=c(1,3), assay=c('7pt_1rep', '10pt_3rep'),
                       sd_prop=0.3)
varying_df <- crossing(sd=1, n=c(50, 200, 800), prop=c(0.05,0.1,0.2), beta=-c(.3, .5, 1), sd_add=c(0.15, 0.5)) %>%
  dplyr::mutate(sim_group=row_number())

#do in parallel with more sims
parallel_sim_df <- crossing(fixed_df, varying_df, sim_rep=c(1:100)) %>%
  dplyr::mutate(sim_unique_id=row_number(),
                batch=sample(1:720, n(), replace = TRUE))
parallel_sim_df

library(batchtools)

bjwork_dir <- file.path(tempdir(), sample(LETTERS, 10) %>% paste(., collapse=''))
bjfiles_dir <- file.path(tempdir(), sample(LETTERS, 10) %>% paste(., collapse=''))
dir.create(bjwork_dir)

#create the registry
reg <- makeRegistry(packages=c('pgxsim', 'tidyverse'), #packages to include
                    work.dir = bjwork_dir,
                    file.dir = bjfiles_dir,
                    make.default=FALSE)

reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = parallel::detectCores())

#map
batchMap(reg=reg, fun=subset_apply, k=1:max(parallel_sim_df$batch),
         more.args=list(df=parallel_sim_df, my_fun=do_simulation_type2, seed=my_seed))
#submitJobs(reg)
submitJobs(reg=reg, resources = list())
getStatus(reg = reg)
done <- findDone(reg=reg)
job_tab <- getJobTable(reg=reg)

#gather results and bind into a dataframe
parallel_res <- reduceResultsList(reg=reg, ids=findDone(reg=reg))
parallel_res_df <- bind_rows(parallel_res) %>%
  inner_join(parallel_sim_df, by='sim_unique_id')

#tidy up
removeRegistry(reg=reg, wait=0)

sess_inf <- sessionInfo()
sess_time <- timestamp(quiet = TRUE)

#save output
save(parallel_res_df, job_tab, fixed_df, varying_df, parallel_sim_df, sess_inf, sess_time, my_seed,
     file = '11_screen_sim.RData')


#estimates
ggplot(parallel_res_df, aes(as.factor(round(beta,2)), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) +
  #facet_grid(prop~sd+n) + ylim(-4,2) +
  facet_grid(sd_add+prop~n) + ylim(-4,2) +
  theme_bw()

#p values
ggplot(parallel_res_df, aes(as.factor(round(beta,2)), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) +
  #facet_grid(prop~sd+n) + ylim(-4,30) +
  facet_grid(sd_add+prop~n) + ylim(-4,20) +
  theme_bw()

#test p values
ggplot(parallel_res_df, aes(as.factor(round(beta,2)), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) +
  #facet_grid(prop~sd+n) + ylim(-4,30) +
  facet_grid(sd_add+prop~n) +
  theme_bw()

#proportion of times that p<0.05
sig_calc <- parallel_res_df %>%
  dplyr::mutate(betap_sig=beta_pval<=0.05,
                testp_sig=test_pval<=0.05,
                rci_sig=abs(1.98*beta_std_err/beta_estimate)<=1) %>%
  dplyr::group_by(sim_group, method, sd, sd_add, n, prop, beta) %>%
  dplyr::summarise(betap_sig=mean(betap_sig),
                   testp_sig=mean(testp_sig),
                   rci_sig=mean(rci_sig))


#power calculations
ggplot(sig_calc, aes(x=n, y=betap_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  facet_grid(prop+sd_add~sd+round(beta,1)) +
  theme_bw()

ggplot(sig_calc, aes(x=n, y=testp_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  facet_grid(prop+sd_add~sd+round(beta,1)) +
  theme_bw()

ggplot(sig_calc, aes(x=n, y=rci_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  facet_grid(prop+sd_add~sd+round(beta,1)) +
  theme_bw()

