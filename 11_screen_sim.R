library(pgxsim)
library(tidyverse)

#set the number of replicate simulations to generate
nreps <- 200
# nreps <- 2 #for testing

#set the seed to make results reproducible
my_seed <- 10001
set.seed(my_seed)

#set output path
outpath <- '11_screen_sim'
if(!dir.exists(outpath)) {
  dir.create(outpath)
}

#define the parameters that will remain fixed
fixed_df <- data_frame(type='discrete', mu=.5, lb=-4, ub=Inf, minconc=0.001, maxconc=30,
                       ndoses=c(7,10), nreps=c(1,3), assay=c('7pt_1rep', '10pt_3rep'),
                       sd_prop=0.3)
varying_df <- crossing(sd=1, n=c(50, 200, 800), prop=c(0.05,0.1,0.2), beta=-c(.3, .5, 1), sd_add=c(0.15, 0.5)) %>%
  dplyr::mutate(sim_group=row_number())

#do in parallel with more sims
nbatches <- max(parallel::detectCores(), nreps*3)  #number of batches to split computation into
parallel_sim_df <- crossing(fixed_df, varying_df, sim_rep=c(1:nreps)) %>%
  dplyr::mutate(sim_unique_id=row_number(),
                batch=sample(1:nbatches, n(), replace = TRUE))
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
     file = file.path(outpath, '11_screen_sim.RData'))


