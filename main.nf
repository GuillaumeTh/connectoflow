#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["apply_t1_labels_transfo":"$params.apply_t1_labels_transfo",
                "output_dir":"$params.output_dir",
                "run_commit":"$params.run_commit",
                "b_thr":"$params.b_thr",
                "nbr_dir":"$params.nbr_dir",
                "ball_stick":"$params.ball_stick",
                "para_diff":"$params.para_diff",
                "perp_diff":"$params.perp_diff",
                "iso_diff":"$params.iso_diff",
                "no_pruning":"$params.no_pruning",
                "no_remove_loops":"$params.no_remove_loops",
                "no_remove_outliers":"$params.no_remove_outliers",
                "min_length":"$params.min_length",
                "max_length":"$params.max_length",
                "loop_max_angle":"$params.loop_max_angle",
                "outlier_threshold":"$params.outlier_threshold",
                "length_weighting":"$params.length_weighting",
                "sh_basis":"$params.sh_basis",
                "run_afd_rd":"$params.run_afd_rd",
                "use_similarity_metric":"$params.use_similarity_metric",
                "nbr_subjects_for_avg_connections":"$params.nbr_subjects_for_avg_connections",
                "processes_register":"$params.processes_register",
                "processes_commit":"$params.processes_commit",
                "processes_afd_rd":"$params.processes_afd_rd",
                "processes_avg_similarity":"$params.processes_avg_similarity",
                "processes_connectivity":"$params.processes_connectivity",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Run Connectivity Construction"
log.info "============================="
log.info ""

log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "[Inputs]"
log.info "Root: $params.input"
log.info "Labels list: $params.labels_list"
log.info ""

root = file(params.input)

Channel.fromPath(file(params.labels_list))
    .into{labels_list_for_compute;labels_list_for_visualize}

metrics_for_compute = Channel
    .fromFilePairs("$root/**/Transform_Metrics/{*_mni.nii.gz}",
                    maxDepth:3,
                    flat: true) {it.parent.parent.name}

metrics_for_compute
    .groupTuple()
    .set{all_metrics_for_compute}

h5_labels_for_compute = Channel
    .fromFilePairs("$root/**/Transform_Data/{*decompose_warped_mni.h5,*labels_warped_mni_int16.nii.gz}",
                    size: 2,
                    maxDepth:3,
                    flat: true) {it.parent.parent.name}

Channel.fromPath("$root/Average_Connections/avg_per_edges")
    .set{edges_for_similarity}

if (params.use_similarity_metric) {
    h5_labels_for_compute.view()
        .join(all_metrics_for_compute)
        .combine(edges_for_similarity)
        .combine(labels_list_for_compute)
        .set{h5_labels_similarity_list_for_compute}
    h5_labels_list_for_compute = Channel.empty()
}
else {
    h5_labels_for_compute
        .join(all_metrics_for_compute)
        .combine(labels_list_for_compute)
        .set{h5_labels_list_for_compute}
    h5_labels_similarity_list_for_compute = Channel.empty()
}

process Compute_Connectivity_with_similiarity {
    cpus params.processes_connectivity
    publishDir = {"${params.output_dir}/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(avg_edges), file(labels_list) from h5_labels_similarity_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_with_similarity

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    for metric in $metrics_list; do
        base_name=\$(basename \${metric/_mni/})
        base_name=\${base_name/_warped/}
        base_name=\${base_name/"${sid}__"/}
        metrics_args="\${metrics_args} --metrics \${metric} \$(basename \$base_name .nii.gz).npy"
    done

    scil_compute_connectivity.py $h5 $labels --force_labels_list $labels_list --volume vol.npy --streamline_count sc.npy \
        --length len.npy --similarity $avg_edges sim.npy \$metrics_args --density_weighting --no_self_connection \
        --include_dps ./ --processes $params.processes_connectivity
    scil_normalize_connectivity.py sc.npy sc_edge_normalized.npy --parcel_volume $labels $labels_list
    scil_normalize_connectivity.py vol.npy sc_vol_normalized.npy --parcel_volume $labels $labels_list
    """
}

process Compute_Connectivity_without_similiarity {
    cpus params.processes_connectivity
    publishDir = {"${params.output_dir}/$sid/Compute_Connectivity"}

    input:
    set sid, file(h5), file(labels), file(metrics), file(labels_list) from h5_labels_list_for_compute

    output:
    set sid, "*.npy" into matrices_for_visualize_without_similarity

    script:
    String metrics_list = metrics.join(", ").replace(',', '')
    """
    metrics_args=""
    for metric in $metrics_list; do
        base_name=\$(basename \${metric/_mni/})
        base_name=\${base_name/_warped/}
        base_name=\${base_name/"${sid}__"/}
        metrics_args="\${metrics_args} --metrics \${metric} \$(basename \$base_name .nii.gz).npy"
    done

    scil_compute_connectivity.py $h5 $labels --force_labels_list $labels_list --volume vol.npy --streamline_count sc.npy \
        --length len.npy \$metrics_args --density_weighting --no_self_connection --include_dps ./ \
        --processes $params.processes_connectivity
    scil_normalize_connectivity.py sc.npy sc_edge_normalized.npy --parcel_volume $labels $labels_list
    scil_normalize_connectivity.py vol.npy sc_vol_normalized.npy --parcel_volume $labels $labels_list
    """
}

matrices_for_visualize_with_similarity
    .concat(matrices_for_visualize_without_similarity)
    .combine(labels_list_for_visualize)
    .set{matrices_labels_list_for_visualize}

process Visualize_Connectivity {
    cpus 1

    input:
    set sid, file(matrices), file(labels_list) from matrices_labels_list_for_visualize

    output:
    set sid, "*.png"

    script:
    String matrices_list = matrices.join(", ").replace(',', '')
    """
    for matrix in $matrices_list; do
        scil_visualize_connectivity.py \$matrix \${matrix/.npy/_matrix.png} --labels_list $labels_list --name_axis \
            --display_legend --histogram \${matrix/.npy/_hist.png} --nb_bins 50 --exclude_zeros --axis_text_size 5 5
    done
    """
}
