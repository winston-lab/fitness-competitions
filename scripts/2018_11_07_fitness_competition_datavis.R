library(tidyverse)
library(broom)
library(magrittr)
library(viridis)
library(ggridges)
library(ggcyto)

label_condition = function(x) paste("A:", x)
label_control = function(x) paste("B:", x)
label_concentration = function(x) paste(x, "mM")
label_days = function(x) paste(x, "days")

theme_default = theme_light() +
    theme(text = element_text(size=8),
          strip.background = element_blank(),
          strip.text = element_text(color="black"),
          strip.text.y = element_text(angle=0, hjust=0))

theme_lineplot = theme_default +
    theme(axis.title.y = element_text(angle=0, hjust=1, vjust=0.5),
          legend.position = c(0.99, 0.07),
          legend.justification = c(1, 0.5),
          panel.grid.minor.x = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.spacing.y = unit(0.2, "cm"),
          legend.key.height=unit(0.5, "cm"),
          legend.title = element_text(margin=margin(0,0,-5,0,"pt")))

main = function(data_path,
                sampling_fraction,
                yfp_cutoff,
                mcherry_cutoff,
                fsc_a_vs_fsc_h_replicate_separate_out = "fsc_a_vs_fsc_h-replicates-separate.png",
                fsc_a_vs_fsc_h_replicate_merged_out = "fsc_a_vs_fsc_h-replicates-merged.png",
                fsc_a_distributions_out = "fsc_a_distributions.png",
                width_distributions_out = "width_distributions.png",
                yfp_distributions_out = "yfp_distributions.png",
                mcherry_distributions_out = "mcherry_distributions.png",
                ratio_lineplot_separate_out = "ratio_lineplot_separate.png",
                proportion_lineplot_separate_out = "proportion_lineplot_separate.png",
                ratio_heatmap_out = "ratio_heatmap.png",
                proportion_heatmap_out = "proportion_heatmap.png",
                ratio_lineplot_merged_out = "ratio_lineplot_merged.png",
                proportion_lineplot_merged_out = "proportion_lineplot_merged.png"){
    df_full = read_tsv(data_path) %>%
        select(treatment_time, diamide, condition, control, condition_fluor, control_fluor,
               replicate, fsc_a=FSC_A, fsc_h=FSC_H, ssc_a=SSC_A, ssc_h=SSC_H,
               yfp=YFP_A, mcherry=mCherry_A, width) %>%
        unite(condition_label, condition, condition_fluor, remove=FALSE) %>%
        unite(control_label, control, control_fluor, remove=FALSE) %>%
        dplyr::filter(fsc_a < max(fsc_a))

    df_sample = df_full %>%
        sample_frac(sampling_fraction)

    # df %>%
    #     filter(fsc_a < max(fsc_a)) %>%
    #     do(lm(.$fsc_h ~ .$fsc_a + I(.$fsc_a^2)) %>%
    #            tidy())
    #
    # m = lm(df$fsc_h ~ df$fsc_a + I(df$fsc_a^2))
    #
    # plot_line = function(x){
    #     m$coefficients[1] +
    #         m$coefficients[2]*x +
    #         m$coefficients[3]*x^2
    # }

    fsc_h_vs_fsc_a = ggplot(data=df_sample,
           aes(x=fsc_a,
               y=fsc_h)) +
        stat_bin_hex(binwidth = c(3e2, 3e2)) +
        # stat_smooth(method="lm", formula=y~x+I(x^2), size=1) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=2),
                           labels=scales::scientific,
                           name="FSC area") +
        scale_y_continuous(breaks=scales::pretty_breaks(n=2),
                           labels=scales::scientific,
                           name="FSC height") +
        scale_fill_viridis(option="inferno",
                           guide=FALSE) +
        theme_default

    ggsave(fsc_a_vs_fsc_h_replicate_separate_out,
           plot = fsc_h_vs_fsc_a +
               facet_grid(control_label + condition_label + treatment_time ~ diamide + replicate,
                          labeller=labeller(control_label=label_control,
                                            condition_label=label_condition,
                                            treatment_time=label_days,
                                            diamide=label_concentration)),
           width=16*3,
           height=9*3,
           units="cm")

    ggsave(fsc_a_vs_fsc_h_replicate_merged_out,
           plot = fsc_h_vs_fsc_a +
               facet_grid(control_label + condition_label + treatment_time ~ diamide,
                          labeller=labeller(control_label=label_control,
                                            condition_label=label_condition,
                                            treatment_time=label_days,
                                            diamide=label_concentration)),
           width=16*2,
           height=9*2,
           units="cm")

    fsc_a_distributions = ggplot(data=df_sample,
           aes(x=fsc_a,
               y=treatment_time,
               group=interaction(replicate, treatment_time))) +
        geom_density_ridges(alpha=0.7,
                            fill="#9ecae1") +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=2),
                           name="FSC area") +
        scale_y_continuous(expand=c(0, 0.03),
                           name="days in diamide") +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller = labeller(condition_label=label_condition,
                                       control_label=label_control,
                                       diamide=label_concentration)) +
        theme_default

    ggsave(fsc_a_distributions_out,
           plot=fsc_a_distributions,
           width=16*2,
           height=9*2,
           units="cm")

    width_distributions = ggplot(data=df_sample,
           aes(x=width,
               y=treatment_time,
               group=interaction(replicate, treatment_time))) +
        geom_density_ridges(alpha=0.7,
                            fill="#9ecae1") +
        scale_x_continuous(expand=c(0,0),
                           breaks=scales::pretty_breaks(n=2),
                           name="pulse width") +
        scale_y_continuous(expand=c(0, 0.03),
                           name="days in diamide") +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller = labeller(condition_label=label_condition,
                                       control_label=label_control,
                                       diamide=label_concentration)) +
        theme_default

    ggsave(width_distributions_out,
           plot=width_distributions,
           width=16*2,
           height=9*2,
           units="cm")

    yfp_distributions = ggplot(data=df_sample,
                               aes(x=yfp,
                                   y=treatment_time,
                                   group=interaction(replicate, treatment_time))) +
        geom_density_ridges(alpha=0.7,
                            fill="#ffeda0") +
        geom_vline(xintercept = yfp_cutoff) +
        scale_x_flowJo_fasinh(expand=c(0,0),
                              name="YFP area",
                              limits=c(NA, quantile(df_sample[["yfp"]], 0.999))) +
        scale_y_continuous(expand=c(0, 0.03),
                           name="days in diamide") +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller = labeller(condition_label=label_condition,
                                       control_label=label_control,
                                       diamide=label_concentration)) +
        theme_default

    ggsave(yfp_distributions_out,
           plot=yfp_distributions,
           width=16*2,
           height=9*2,
           units="cm")

    mcherry_distributions = ggplot(data=df_sample,
                                   aes(x=mcherry,
                                       y=treatment_time,
                                       group=interaction(replicate, treatment_time))) +
        geom_density_ridges(alpha=0.7,
                            fill="#ef3b2c") +
        geom_vline(xintercept = mcherry_cutoff) +
        scale_x_flowJo_fasinh(expand=c(0,0),
                              name="mCherry area",
                              limits=c(NA, quantile(df_sample[["mcherry"]], 0.999))) +
        scale_y_continuous(expand=c(0, 0.03),
                           name="days in diamide") +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller = labeller(condition_label=label_condition,
                                       control_label=label_control,
                                       diamide=label_concentration)) +
        theme_default

    ggsave(mcherry_distributions_out,
           plot=mcherry_distributions,
           width=16*2,
           height=9*2,
           units="cm")

    df_full %<>%
        mutate(fluor_identity=case_when(mcherry > mcherry_cutoff & yfp <= yfp_cutoff ~ "mcherry",
                                        mcherry <= mcherry_cutoff & yfp > yfp_cutoff ~ "yfp",
                                        mcherry > mcherry_cutoff & yfp > yfp_cutoff ~ "doublet",
                                        mcherry <= mcherry_cutoff & yfp <= yfp_cutoff ~ "non-fluorescing"))


    df_summary = df_full %>%
        dplyr::filter(fluor_identity %in% c("mcherry", "yfp")) %>%
        group_by(treatment_time, diamide, replicate,
                 condition_label, condition, condition_fluor,
                 control_label, control, control_fluor) %>%
        count(fluor_identity) %>%
        spread(key=fluor_identity, value=n) %>%
        mutate(ratio=if_else(condition_fluor=="mCherry",
                             log2(mcherry/yfp),
                             log2(yfp/mcherry)),
               proportion=if_else(condition_fluor=="mCherry",
                                  mcherry/(mcherry+yfp),
                                  yfp/(mcherry+yfp)),
               experiment_type=if_else(condition==control,
                                       "self vs. self",
                                       "mutant vs. WT"))

    ratio_heatmap = ggplot(data=df_summary) +
        geom_raster(aes(x=replicate, y=treatment_time, fill=ratio)) +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller=labeller(control_label = label_control,
                                     condition_label = label_condition,
                                     diamide = label_concentration)) +
        scale_x_continuous(expand=c(0,0),
                           name="replicate",
                           breaks=1:3) +
        scale_y_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_fill_distiller(palette="PRGn",
                             limits=c(-3,3),
                             labels=c("≤-3", -2, -1, 0, 1, 2, "≥3"),
                             name=expression(log[2] ~ frac(A, B)),
                             oob=scales::squish,
                             guide=guide_colorbar(barwidth=unit(5, "cm"),
                                                  barheight=unit(0.5, "cm"),
                                                  title.vjust=0.75)) +
        theme_default +
        theme(panel.border=element_blank(),
              panel.grid=element_blank(),
              legend.position="top")
    ggsave(ratio_heatmap_out,
           plot=ratio_heatmap,
           width=16*1.5,
           height=9*1.5,
           units="cm")

    proportion_heatmap = ggplot(data=df_summary) +
        geom_raster(aes(x=replicate, y=treatment_time, fill=proportion)) +
        facet_grid(control_label + condition_label ~ diamide,
                   labeller=labeller(control_label = label_control,
                                     condition_label = label_condition,
                                     diamide = label_concentration)) +
        scale_x_continuous(expand=c(0,0),
                           name="replicate",
                           breaks=1:3) +
        scale_y_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_fill_distiller(palette="PRGn",
                             limits=c(0,1),
                             name="proportion \'A\' cells",
                             oob=scales::squish,
                             guide=guide_colorbar(barwidth=unit(5, "cm"),
                                                  barheight=unit(0.5, "cm"),
                                                  title.vjust=0.75)) +
        theme_default +
        theme(panel.border=element_blank(),
              panel.grid=element_blank(),
              legend.position="top")

    ggsave(proportion_heatmap_out,
           plot=proportion_heatmap,
           width=16*1.5,
           height=9*1.5,
           units="cm")

    ratio_lineplot_separate = ggplot(data=df_summary) +
        geom_hline(yintercept=0) +
        geom_line(aes(x=treatment_time,
                      y=ratio,
                      linetype=experiment_type,
                      color=condition_fluor,
                      group=interaction(replicate, condition_label, control_label)),
                  alpha=0.6) +
        scale_x_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_y_continuous(name=expression(log[2] ~ frac(A, B)),
                           breaks=scales::pretty_breaks(n=3)) +
        scale_linetype_manual(values=c("solid", "dotted"),
                              name="A vs. B") +
        scale_color_brewer(palette = "Set1",
                           name="\'A\' fluorophore") +
        facet_wrap(~diamide, ncol=4,
                   labeller = labeller(diamide = label_concentration)) +
        theme_lineplot

    ggsave(ratio_lineplot_separate_out,
           plot=ratio_lineplot_separate,
           width=16, height=9, units="cm")

    proportion_lineplot_separate = ggplot(data=df_summary) +
        geom_hline(yintercept=0.5) +
        geom_line(aes(x=treatment_time,
                      y=proportion,
                      linetype=experiment_type,
                      color=condition_fluor,
                      group=interaction(replicate, condition_label, control_label)),
                  alpha=0.6) +
        scale_x_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_y_continuous(name="proportion\n\'A\' cells",
                           breaks=scales::pretty_breaks(n=3)) +
        scale_linetype_manual(values=c("solid", "dotted"),
                              name="A vs. B") +
        scale_color_brewer(palette = "Set1",
                           name="\'A\' fluorophore") +
        facet_wrap(~diamide, ncol=4,
                   labeller = labeller(diamide = label_concentration)) +
        theme_lineplot

    ggsave(proportion_lineplot_separate_out,
           plot=proportion_lineplot_separate,
           width=16, height=9, units="cm")

    df_summary_merged = df_summary %>%
        group_by(treatment_time, diamide, experiment_type) %>%
        summarise(ratio_mean=mean(ratio),
                  ratio_sd=sd(ratio),
                  proportion_mean=mean(proportion),
                  proportion_sd=sd(proportion))

    ratio_lineplot_merged = ggplot(data=df_summary_merged) +
        geom_hline(yintercept=0) +
        geom_ribbon(aes(x=treatment_time,
                        ymax=ratio_mean+ratio_sd,
                        ymin=ratio_mean-ratio_sd,
                        group=experiment_type,
                        alpha=experiment_type),
                    fill="#114477") +
        geom_line(aes(x=treatment_time,
                      y=ratio_mean,
                      linetype=experiment_type,
                      alpha=experiment_type),
                  color="#114477") +
        scale_x_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_y_continuous(name=expression(log[2] ~ frac(A, B)),
                           breaks=scales::pretty_breaks(n=3)) +
        scale_alpha_manual(values=c(0.6, 0.3),
                           name="A vs. B") +
        scale_linetype_manual(values=c("solid", "dotted"),
                              name="A vs. B") +
        facet_wrap(~diamide, ncol=4,
                   labeller = labeller(diamide = label_concentration)) +
        theme_lineplot

    ggsave(ratio_lineplot_merged_out,
           plot=ratio_lineplot_merged,
           width=16, height=9, units="cm")

    proportion_lineplot_merged = ggplot(data=df_summary_merged) +
        geom_hline(yintercept=0.5) +
        geom_ribbon(aes(x=treatment_time,
                        ymax=proportion_mean+proportion_sd,
                        ymin=proportion_mean-proportion_sd,
                        group=experiment_type,
                        alpha=experiment_type),
                    fill="#114477") +
        geom_line(aes(x=treatment_time,
                      y=proportion_mean,
                      linetype=experiment_type,
                      alpha=experiment_type),
                  color="#114477") +
        scale_x_continuous(expand=c(0,0),
                           name="days in diamide",
                           breaks=0:2) +
        scale_y_continuous(name="proportion\n\'A\' cells",
                           breaks=scales::pretty_breaks(n=3)) +
        scale_alpha_manual(values=c(0.6, 0.3),
                           name="A vs. B") +
        scale_linetype_manual(values=c("solid", "dotted"),
                              name="A vs. B") +
        facet_wrap(~diamide, ncol=4,
                   labeller = labeller(diamide = label_concentration)) +
        theme_lineplot

    ggsave(proportion_lineplot_merged_out,
           plot=proportion_lineplot_merged,
           width=16, height=9, units="cm")
}

main(data_path=snakemake@input[["data_path"]],
     sampling_fraction=snakemake@params[["sampling_fraction"]],
     yfp_cutoff=snakemake@params[["yfp_cutoff"]],
     mcherry_cutoff=snakemake@params[["mcherry_cutoff"]],
     fsc_a_vs_fsc_h_replicate_separate_out=snakemake@output[["fsc_a_vs_fsc_h_replicate_separate_out"]],
     fsc_a_vs_fsc_h_replicate_merged_out=snakemake@output[["fsc_a_vs_fsc_h_replicate_merged_out"]],
     fsc_a_distributions_out=snakemake@output[["fsc_a_distributions_out"]] ,
     width_distributions_out=snakemake@output[["width_distributions_out"]] ,
     yfp_distributions_out=snakemake@output[["yfp_distributions_out"]],
     mcherry_distributions_out=snakemake@output[["mcherry_distributions_out"]] ,
     ratio_lineplot_separate_out=snakemake@output[["ratio_lineplot_separate_out"]],
     proportion_lineplot_separate_out=snakemake@output[["proportion_lineplot_separate_out"]],
     ratio_heatmap_out=snakemake@output[["ratio_heatmap_out"]],
     proportion_heatmap_out=snakemake@output[["proportion_heatmap_out"]],
     ratio_lineplot_merged_out=snakemake@output[["ratio_lineplot_merged_out"]],
     proportion_lineplot_merged_out=snakemake@output[["proportion_lineplot_merged_out"]])

