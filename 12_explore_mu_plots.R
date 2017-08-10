library(pgxsim)
library(tidyverse)

#define the output path and load the simulation data (if present)
outpath <- '12_explore_mu'
res_file <- file.path(outpath, '12_explore_mu.RData')

if(file.exists(res_file)) {
  load(res_file)
} else {
  stop('Need to run 12_explore_mu.R script first')
}

#get rid of truncated nls_lm
parallel_res_df <- dplyr::filter(parallel_res_df, method!='nls_lm_t0')

#show simulation setup
DT::datatable(varying_df)
DT::datatable(fixed_df)

#beta estimates 7pt_1rep
beta_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
ggplot(aes(as.factor(mu), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 3A: Explore mu 7pt_1rep betas') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both) +
  ylim(-2,4) + xlab('mu') +
  theme_bw()
beta_7pt_1rep
ggsave('12_beta_7pt_1rep_fig3a.png', beta_7pt_1rep, path=outpath, width = 8, height=6)


#beta estimates 10pt_3rep
beta_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(mu), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 3B: Explore mu 10pt_3rep betas') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both) +
  ylim(-2,4) + xlab('mu') +
  theme_bw()
beta_10pt_3rep
ggsave('12_beta_10pt_3rep_fig3b.png', beta_10pt_3rep, path=outpath, width = 8, height=6)


#beta pval 7pt_1rep
betapval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(mu), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Explore mu 7pt_1rep beta pval') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('mu')
betapval_7pt_1rep
ggsave('12_betapval_7pt_1rep.png', betapval_7pt_1rep, path=outpath, width = 8, height=6)


#beta pval 10pt_3rep
betapval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(mu), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Explore mu 10pt_3rep beta pval') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('mu')
betapval_10pt_3rep
ggsave('12_betapval_10pt_3rep.png', betapval_10pt_3rep, path=outpath, width = 8, height=6)


#test pval 7pt_1rep
testpval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(mu), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 3C: Explore mu 7pt_1rep test pval') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('mu')
testpval_7pt_1rep
ggsave('12_testpval_7pt_1rep_fig3c.png', testpval_7pt_1rep, path=outpath, width = 8, height=6)


#test pval 10pt_3rep
testpval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(mu), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 3D: Explore mu 10pt_3rep test pval') +
  facet_grid(assay+sd~n+prop+sd_add, labeller = label_both, scales='free_y') +
  theme_bw() + xlab('mu')
testpval_10pt_3rep
ggsave('12_testpval_10pt_3rep_fig3d.png', testpval_10pt_3rep, path=outpath, width = 8, height=6)

sess_time
my_seed
sess_inf