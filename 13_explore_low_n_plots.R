library(pgxsim)
library(tidyverse)


setwd('~/Documents/Work/2017/20170531 Dose Response Simulation/aws_runs/')
load('13_explore_low_n.RData')

#get rid of truncated nls_lm
parallel_res_df <- dplyr::filter(parallel_res_df, method!='nls_lm_t0')

#show simulation setup
DT::datatable(varying_df)
DT::datatable(fixed_df)

#beta estimates 7pt_1rep
beta_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(beta), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 4A: Explore low n 7pt_1rep betas') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw() + ylim(-4,2) + xlab('beta')
beta_7pt_1rep
ggsave('13_beta_7pt_1rep_fig4a.png', beta_7pt_1rep, path='13_explore_low_n_plots/', width = 10, height=6)


#beta estimates 10pt_3rep
beta_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(beta), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 4B: Explore low n 10pt_3rep betas') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw() + ylim(-4,2) + xlab('beta')
beta_10pt_3rep
ggsave('13_beta_10pt_3rep_fig4b.png', beta_10pt_3rep, path='13_explore_low_n_plots/', width = 10, height=6)


#beta pvals 7pt_1rep
betapval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(beta), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Explore low n 7pt_1rep beta pvals') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw()
betapval_7pt_1rep

#beta pval 10pt_3rep
betapval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(beta), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Explore low n 10pt_3rep beta pvals') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw()
betapval_10pt_3rep

#test pvals 7pt_1rep
testpval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(beta), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 4C: Explore low n 7pt_1rep test pvals') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('beta')
testpval_7pt_1rep
ggsave('13_testpval_7pt_1rep_fig4c.png', testpval_7pt_1rep, path='13_explore_low_n_plots/', width = 10, height=6)


#test pval 10pt_3rep
testpval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(beta), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 4D: Explore low n 10pt_3rep test pvals') +
  facet_grid(assay+sd_add~sd+n+prop, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('beta')
testpval_10pt_3rep
ggsave('13_testpval_10pt_3rep_fig4d.png', testpval_10pt_3rep, path='13_explore_low_n_plots/', width = 10, height=6)

sess_time
my_seed
sess_inf

