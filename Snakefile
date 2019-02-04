#!/usr/bin/python

configfile: "config.yaml"

EXPERIMENTS = config["experiments"]

rule target:
    input:
        "config.yaml",
        expand("parsed_data/{experiment}_fitness_competition_flow_cytometry.tsv.gz", experiment=EXPERIMENTS),
        expand("{experiment}_ratio_heatmap.svg", experiment=EXPERIMENTS)

rule parse_flow_cytometry_data:
    input:
        samplesheet = lambda wc: EXPERIMENTS[wc.experiment]["samplesheet"]
    output:
        "parsed_data/{experiment}_fitness_competition_flow_cytometry.tsv.gz"
    params:
        uncompressed = "parsed_data/{experiment}_fitness_competition_flow_cytometry.tsv"
    conda:
        "envs/parse_flow_cytometry.yaml"
    shell: """
        python scripts/parse_flow_cytometry.py -i {input.samplesheet} -o {params.uncompressed}
        pigz -f {params.uncompressed} > {output}
        """

rule fitness_competition_datavis:
    input:
        data_path = "parsed_data/{experiment}_fitness_competition_flow_cytometry.tsv.gz"
    output:
        fsc_a_vs_fsc_h_replicate_separate_out = "{experiment}_fsc_a_vs_fsc_h-replicates-separate.svg",
        fsc_a_vs_fsc_h_replicate_merged_out = "{experiment}_fsc_a_vs_fsc_h-replicates-merged.svg",
        fsc_a_distributions_out = "{experiment}_fsc_a_distributions.svg",
        width_distributions_out = "{experiment}_width_distributions.svg",
        yfp_distributions_out = "{experiment}_yfp_distributions.svg",
        mcherry_distributions_out = "{experiment}_mcherry_distributions.svg",
        ratio_lineplot_separate_out = "{experiment}_ratio_lineplot_separate.svg",
        proportion_lineplot_separate_out = "{experiment}_proportion_lineplot_separate.svg",
        ratio_heatmap_out = "{experiment}_ratio_heatmap.svg",
        proportion_heatmap_out = "{experiment}_proportion_heatmap.svg",
        ratio_lineplot_merged_out = "{experiment}_ratio_lineplot_merged.svg",
        proportion_lineplot_merged_out = "{experiment}_proportion_lineplot_merged.svg"
    params:
        sampling_fraction = lambda wc: EXPERIMENTS[wc.experiment]["sampling_fraction"],
        yfp_cutoff = lambda wc: EXPERIMENTS[wc.experiment]["yfp_cutoff"],
        mcherry_cutoff = lambda wc: EXPERIMENTS[wc.experiment]["mcherry_cutoff"],
    conda:
        "envs/fitness_competition_datavis.yaml"
    script:
        "scripts/2018_11_07_fitness_competition_datavis.R"




