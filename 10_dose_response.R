library(tidyverse)
library(pgxsim)

set.seed(10001)

#set output path
outpath <- '10_dose_response'
if(!dir.exists(outpath)) {
  dir.create(outpath)
}

#simulate cell lines and dose response data
pIC50_seq <- seq(-5,4,0.05)
pIC50_data <- data_frame(cell_id=1:length(pIC50_seq), gene=0, pIC50=pIC50_seq)
param_df <- tidyr::crossing(minconc=0.001, maxconc=30, ndoses=10,
                            nreps=3, sd_prop=c(0.1, 0.2, 0.4, 0.8),
                            sd_add=c(0.1, 0.2, 0.4, 0.8)) %>%
  dplyr::mutate(sim_no=row_number())

fixed_df <- data_frame(type='discrete', minconc=0.001, maxconc=30,
                       ndoses=c(7,10), nreps=c(1,3), assay=c('7pt_1rep', '10pt_3rep'))
varying_df <- crossing(sd_prop=c(0.1, 0.2, 0.4, 0.8),
                       sd_add=c(0.1, 0.2, 0.4, 0.8))
param_df <- tidyr::crossing(fixed_df, varying_df) %>%
  dplyr::mutate(sim_no=row_number())

simulation_df <- tidyr::crossing(pIC50_data, param_df)
simulation_df
simulation_df <-  simulation_df %>% dplyr::mutate(
  dr_data=purrr::pmap(.l=list(pIC50=pIC50, minconc=minconc, maxconc=maxconc, ndoses=ndoses,
                              nreps=nreps, sd_prop=sd_prop, sd_add=sd_add),
                      .f=sim_dose_response)) %>%
  dplyr::select(sim_no, cell_id, assay, gene, pIC50, sd_prop, sd_add, dr_data)
simulation_df %>% dplyr::select(sim_no, cell_id, dr_data)

#show simulated dose response data for cell line 100
cl100_dr_plot <- simulation_df %>%
  dplyr::filter(cell_id==100 & assay=='10pt_3rep') %>%
  tidyr::unnest() %>%
  ggplot(aes(conc,resp)) + geom_point() + scale_x_log10() + ylim(-1,2) +
  theme_bw() + facet_wrap(sd_prop~sd_add, labeller = label_both) + ggtitle('Fig 1: Simulated Dose Response Curves') 
cl100_dr_plot
ggsave('cl100_dr_plot_fig1.png', cl100_dr_plot, path=outpath, width = 8, height=6)

#fit dose response curves using nls
nls_fits <- simulation_df %>%
  dplyr::mutate(fit=purrr::map(dr_data, nls_fit))

#extract the pIC50s but truncate at max/min pIC50 +/- 3
nls_results <- nls_fits %>%
  dplyr::mutate(res=purrr::map(fit, nls_extract, upper_trunc=3, lower_trunc=3)) %>%
  dplyr::select(-dr_data, -fit) %>%
  tidyr::unnest()

nls_fit_plot_10pt <- nls_results %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(pIC50, nls_pIC50)) +
  geom_point(shape=21) + ggtitle('Fig 2C: NLS 10pt_3rep pIC50s') +
  geom_abline(slope=1, intercept=0, colour='blue', linetype='dotted') +
  geom_vline(xintercept=log10(c(0.001, 30)), colour='darkgreen', linetype='dashed') +
  geom_smooth(method='lm', colour='red') +
  facet_grid(assay+sd_prop~sd_add, labeller = label_both) +
  theme_bw() + xlim(-6,5) + ylim(-6,5)
nls_fit_plot_10pt
ggsave('nls_fit_plot_10pt_fig2c.png', nls_fit_plot_10pt, path=outpath, width = 8, height=6)

nls_fit_plot_7pt <- nls_results %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(pIC50, nls_pIC50)) +
  geom_point(shape=21) + ggtitle('Fig 2A: NLS 7pt_1rep pIC50s') +
  geom_abline(slope=1, intercept=0, colour='blue', linetype='dotted') +
  geom_vline(xintercept=log10(c(0.001, 30)), colour='darkgreen', linetype='dashed') +
  geom_smooth(method='lm', colour='red') +
  facet_grid(assay+sd_prop~sd_add, labeller = label_both) +
  theme_bw() + xlim(-6,5) + ylim(-6,5)
nls_fit_plot_7pt
ggsave('nls_fit_plot_7pt_fig2a.png', nls_fit_plot_7pt, path=outpath, width = 8, height=6)


#fit dose response curves using nlme
nlme_results <- simulation_df %>%
  tidyr::unnest() %>%
  dplyr::group_by(sim_no) %>%
  tidyr::nest() %>%
  dplyr::mutate(fit=purrr::map(data, nlme_fit),
                res=purrr::map(fit, nlme_extract)) %>%
  dplyr::select(-data, -fit) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(simulation_df, by=c('sim_no', 'cell_id'))

nlme_fit_plot_10pt <- nlme_results %>%
  dplyr::filter(assay=='10pt_3rep') %>%
  ggplot(aes(pIC50, nlme_pIC50)) +
  geom_point(shape=21) +ggtitle('Fig 2D: NLME 10pt_3rep pIC50s') +
  geom_abline(slope=1, intercept=0, colour='blue', linetype='dotted') +
  geom_vline(xintercept=log10(c(0.001, 30)), colour='darkgreen', linetype='dashed') +
  geom_smooth(method='lm', colour='red') +
  facet_grid(assay+sd_prop~sd_add, labeller = label_both) +
  theme_bw() + xlim(-6,5) + ylim(-6,5)
nlme_fit_plot_10pt
ggsave('nlme_fit_plot_10pt_fig2d.png', nlme_fit_plot_10pt, path=outpath, width = 8, height=6)

nlme_fit_plot_7pt <- nlme_results %>%
  dplyr::filter(assay=='7pt_1rep') %>%
  ggplot(aes(pIC50, nlme_pIC50)) +
  geom_point(shape=21) + ggtitle('Fig 2B: NLME 7pt_1rep pIC50s') +
  geom_abline(slope=1, intercept=0, colour='blue', linetype='dotted') +
  geom_vline(xintercept=log10(c(0.001, 30)), colour='darkgreen', linetype='dashed') +
  geom_smooth(method='lm', colour='red') +
  facet_grid(assay+sd_prop~sd_add, labeller = label_both) +
  theme_bw() + xlim(-6,5) + ylim(-6,5)
nlme_fit_plot_7pt
ggsave('nlme_fit_plot_7pt_fig2b.png', nlme_fit_plot_7pt, path=outpath, width = 8, height=6)
