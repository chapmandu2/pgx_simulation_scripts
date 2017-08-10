library(pgxsim)
library(tidyverse)

#define the output path and load the simulation data (if present)
outpath <- '11_screen_sim'
res_file <- file.path(outpath, '11_screen_sim.RData')

if(file.exists(res_file)) {
  load(res_file)
} else {
  stop('Need to run 11_screen_sim.R script first')
}

#get rid of truncated nls_lm
parallel_res_df <- dplyr::filter(parallel_res_df, method!='nls_lm_t0')

#show simulation setup
DT::datatable(varying_df)
DT::datatable(fixed_df)

#fix test p-value for nlme_gene
parallel_res_df <- parallel_res_df %>%
  dplyr::mutate(test_pval = case_when(
    method == 'nlme_gene' ~ stats::pchisq(test_stat,1,lower.tail=FALSE),
    TRUE ~ test_pval
  ))

#beta estimates 7pt_1rep
beta_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
ggplot(aes(as.factor(n), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 5A: Screen sim 7pt_1rep betas') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + ylim(-3,2) + xlab('Number of cell lines')
beta_7pt_1rep
ggsave('11_beta_7pt_1rep_fig5a.png', beta_7pt_1rep, path=outpath, width = 10, height=6)


#beta estimates 10pt_3rep
beta_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(n), beta_estimate, colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 5B: Screen sim 10pt_3rep betas') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + ylim(-3,2) + xlab('Number of cell lines')
beta_10pt_3rep
ggsave('11_beta_10pt_3rep_fig5b.png', beta_10pt_3rep, path=outpath, width = 10, height=6)


#beta pval 7pt_1rep
betapval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(n), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Screen sim 7pt_1rep beta pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw()  + xlab('Number of cell lines')
betapval_7pt_1rep

#beta pval 10pt_3rep
betapval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(n), -log10(beta_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Screen sim 10pt_3rep beta pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw()  + xlab('Number of cell lines')
betapval_10pt_3rep

#test pval 7pt_1rep
testpval_7pt_1rep <- parallel_res_df %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(as.factor(n), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 5C: Screen sim 7pt_1rep test pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw()  + xlab('Number of cell lines')
testpval_7pt_1rep
ggsave('11_testpval_7pt_1rep_fig5c.png', testpval_7pt_1rep, path=outpath, width = 10, height=6)


#test pval 10pt_3rep
testpval_10pt_3rep <- parallel_res_df %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(as.factor(n), -log10(test_pval), colour=method)) +
  geom_boxplot(outlier.size = 0) + ggtitle('Fig 5D: Screen sim 10pt_3rep test pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + xlab('Number of cell lines')
testpval_10pt_3rep
ggsave('11_testpval_10pt_3rep_fig5d.png', testpval_10pt_3rep, path=outpath, width = 10, height=6)


# Power calculation -------------------------------------------------------


#proportion of times that p<0.05
sig_calc <- parallel_res_df %>%
  dplyr::mutate(betap_sig=beta_pval<=0.05,
                testp_sig=test_pval<=0.05,
                rci_sig=abs(1.98*beta_std_err/beta_estimate)<=1) %>%
  dplyr::group_by(sim_group, assay, method, sd, sd_add, n, prop, beta) %>%
  dplyr::summarise(betap_sig=mean(betap_sig),
                   testp_sig=mean(testp_sig),
                   rci_sig=mean(rci_sig))


#power calc plots
#beta pval power 7pt_1rep
pwr_betapval_7pt_1rep <-sig_calc %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(x=n, y=betap_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Screen sim 7pt_1rep beta pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + xlab('Number of cell lines')
pwr_betapval_7pt_1rep

#beta pval power 10pt_3rep
pwr_betapval_10pt_3rep <-sig_calc %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(x=n, y=betap_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Screen sim 10pt_3rep beta pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + xlab('Number of cell lines')
pwr_betapval_10pt_3rep


#test pval power 7pt_1rep
pwr_testpval_7pt_1rep <-sig_calc %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(x=n, y=testp_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Fig 6A: Screen sim 7pt_1rep test pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + xlab('Number of cell lines')
pwr_testpval_7pt_1rep
ggsave('11_pwr_testpval_7pt_1rep_fig6a.png', pwr_testpval_7pt_1rep, path=outpath, width = 10, height=6)


#test pval power 10pt_3rep
pwr_testpval_10pt_3rep <-sig_calc %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(x=n, y=testp_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Fig 6B: Screen sim 10pt_3rep test pval') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() + xlab('Number of cell lines')
pwr_testpval_10pt_3rep
ggsave('11_pwr_testpval_10pt_3rep_fig6b.png', pwr_testpval_10pt_3rep, path=outpath, width = 10, height=6)


#rci power 7pt_1rep
pwr_rci_7pt_1rep <-sig_calc %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(x=n, y=testp_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Screen sim 7pt_1rep rci') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw()
pwr_rci_7pt_1rep

#rci power 10pt_3rep
pwr_rci_10pt_3rep <-sig_calc %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(x=n, y=testp_sig, colour=method, linetype=method)) +
  geom_line(alpha=0.8, size=1) +
  ggtitle('Screen sim 10pt_3rep rci') +
  facet_grid(assay+prop~sd+sd_add+beta, labeller = label_both) +
  theme_bw() 
pwr_rci_10pt_3rep

sess_time
my_seed
sess_inf
