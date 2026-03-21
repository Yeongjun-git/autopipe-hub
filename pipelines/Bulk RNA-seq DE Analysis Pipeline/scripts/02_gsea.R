# 중복 gene 제거 유틸리티: 같은 gene_id가 여러 번 있으면 절대값이 큰 것을 선택
deduplicate_named_vector <- function(genes) {
  df <- data.frame(name = names(genes), value = genes, stringsAsFactors = FALSE)
  df <- df[!is.na(df$name) & df$name != "", ]
  df$abs_value <- abs(df$value)
  df <- df[order(-df$abs_value), ]
  df <- df[!duplicated(df$name), ]
  result <- setNames(df$value, df$name)
  sort(result, decreasing = TRUE)
}

convert_vector_to_entrez <- function(genes, orgdb, keytype) {
  original_ids <- names(genes)
  conversion <- bitr(original_ids, fromType = keytype, toType = "ENTREZID", OrgDb = orgdb)
  matched <- conversion[conversion[[keytype]] %in% original_ids, ]
  matched$log2FC <- genes[matched[[keytype]]]
  entrez_vector <- setNames(matched$log2FC, matched$ENTREZID)
  entrez_vector <- deduplicate_named_vector(entrez_vector)
  return(entrez_vector)
}

empty_gsea_df <- function() {
  data.frame(
    ID = character(), Description = character(), setSize = character(),
    enrichmentScore = character(), NES = character(), pvalue = numeric(),
    p.adjust = numeric(), qvalue = numeric(), rank = character(),
    leading_edge = character(), core_enrichment = character(), stringsAsFactors = FALSE
  )
}

description_gsea_df <- function() {
  data.frame(
    Column = c("Reference","ID","Description","setSize","enrichmentScore","NES","pvalue","p.adjust","qvalue","rank","leading_edge","core_enrichment"),
    Meaning = c(
      "https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html",
      "Term ID","Term name","Number of genes in the term",
      "Raw enrichment score (ES)","Normalized Enrichment Score",
      "The p-value of the ES is calculated using a permutation test",
      "Adjusted p-value (Benjamini-Hochberg method)",
      "Adjusted p-value (Estimated false discovery rate (FDR))",
      "Rank in the ordered gene list where the enrichment score peak occurs",
      "Summary of the leading edge subset",
      "The actual genes contributing to the enrichment signal"
    )
  )
}

filter_gsea_result <- function(df, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
  df %>% arrange(p.adjust) %>% filter(p.adjust < pvalue_cutoff, setSize >= min_set_size, setSize <= max_set_size)
}

adjust_plot_size <- function(top_df) {
  width <- 8; height <- 8
  max_label_length <- max(nchar(as.character(top_df$Description)))
  counts = nrow(top_df)
  list(width = width + max_label_length * 0.1, height = height - (20 - counts) * 0.3)
}

select_top_GSEA_abs <- function(GSEA_result, top_n = 30) {
  Len = nrow(GSEA_result)
  GSEA_result %>% data.frame() %>% arrange(p.adjust) %>% head(n = min(top_n, Len)) %>% arrange(NES) -> top_pathway_data
  top_pathway_data$description = paste0(top_pathway_data$Description, " (setsize: ", top_pathway_data$setSize, ") ")
  factor(top_pathway_data$description, levels = top_pathway_data$description) -> top_pathway_data$description
  return(top_pathway_data)
}

draw_gsea_summary_plot <- function(GSEA_result, prefix_suffix, path, Comparison_name, top_n = 30, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
  set.seed(116)
  sort_df <- GSEA_result %>% arrange(p.adjust)
  sort_df = filter_gsea_result(sort_df, pvalue_cutoff, min_set_size, max_set_size)
  top_df = select_top_GSEA_abs(sort_df, top_n)
  plot_size <- adjust_plot_size(top_df)
  pdf(file = paste0(path, "3.GSEA/1.figure/", Comparison_name, "_", prefix_suffix, ".pdf"), width = plot_size$width, height = plot_size$height)
  p = ggplot(top_df) +
    geom_bar(aes(x = description, y = NES, fill = p.adjust), stat = 'identity') +
    scale_fill_continuous(low = "red", high = "blue") +
    coord_flip(ylim = c(-3.2, 3.2)) + scale_y_continuous(breaks = seq(-3, 3, 1)) +
    theme_bw(14) + labs(title = paste0(Comparison_name, "_", prefix_suffix), x = "GO Term") +
    theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())
  print(p); dev.off()
}

perform_go_gsea <- function(genes, keytype, orgdb, Comparison_name, path) {
  ontologies <- c("BP", "CC", "MF"); go_results <- list()
  for (ont in ontologies) {
    ego <- gseGO(gene = genes, OrgDb = orgdb, ont = ont, keyType = keytype, pAdjustMethod = "BH")
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      result_sorted <- as.data.frame(ego)[order(as.data.frame(ego)$p.adjust), ]
    } else { result_sorted <- empty_gsea_df() }
    go_results[[ont]] <- result_sorted
  }
  return(go_results)
}

perform_kegg_gsea <- function(genes, keytype, orgdb, Org_code) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)
  ekegg <- gseKEGG(gene = conv, organism = Org_code, keyType = "ncbi-geneid", pAdjustMethod = "BH")
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)[order(as.data.frame(ekegg_readable)$p.adjust), ]
  } else { result <- empty_gsea_df() }
  return(result)
}

perform_reactome_gsea <- function(genes, keytype, orgdb, organism) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)
  erct <- gsePathway(gene = conv, organism = organism, pAdjustMethod = "BH")
  if (!is.null(erct) && nrow(as.data.frame(erct)) > 0) {
    readable <- setReadable(erct, OrgDb = orgdb, keyType = "ENTREZID")
    result_sorted <- as.data.frame(readable)[order(as.data.frame(readable)$p.adjust), ]
  } else { result_sorted <- empty_gsea_df() }
  return(result_sorted)
}

perform_custom_geneset_gsea <- function(genes, Custom_geneset_gmt) {
  if (length(genes) == 0) return(empty_gsea_df())
  Custom_df <- read.gmt(Custom_geneset_gmt)
  gsea_res <- GSEA(geneList = genes, TERM2GENE = Custom_df, pAdjustMethod = "BH")
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    return(as.data.frame(gsea_res)[order(as.data.frame(gsea_res)$p.adjust), ])
  } else { return(empty_gsea_df()) }
}

run_specific_gsea_analysis <- function(path, data, Comparison_name, specific_geneset, pvalue_cutoff = 1) {
  try({
    if (is.null(specific_geneset)) return(invisible(NULL))
    Genes <- setNames(data$log2FoldChange, data$gene_id)
    Genes <- deduplicate_named_vector(Genes)
    hallmark_df <- read.gmt(specific_geneset)
    gsea_res <- GSEA(geneList = Genes, TERM2GENE = hallmark_df, pAdjustMethod = "BH", pvalueCutoff = pvalue_cutoff)
    draw_Specific_GSEA_plot(path = path, data = gsea_res, Name = Comparison_name)
  })
}

GSEA_analysis = function(path, data, Organism, Comparison_name, Keytype, Top_n = 20, Custom_geneset_gmt = NULL, prefix_suffix = "GO_BP") {
  make_folder(path, "3.GSEA/"); make_folder(path, "3.GSEA/0.excel/"); make_folder(path, "3.GSEA/1.figure/")
  if (Organism == "human") {
    OrgDb = org.Hs.eg.db; Org_code = "hsa"
  } else if (Organism == "mouse") {
    OrgDb = org.Mm.eg.db; Org_code = "mmu"
  }
  Hallmark <- resolve_hallmark_path(Organism)

  all_results = list()
  Genes <- setNames(data$log2FoldChange, data$gene_id)
  Genes <- deduplicate_named_vector(Genes)
  prefix <- Comparison_name
  go_res <- perform_go_gsea(Genes, Keytype, OrgDb)
  kegg_res <- perform_kegg_gsea(Genes, Keytype, OrgDb, Org_code)
  reactome_res <- perform_reactome_gsea(Genes, Keytype, OrgDb, Organism)
  hallmark_res <- if (Hallmark != "" && file.exists(Hallmark)) {
    perform_custom_geneset_gsea(Genes, Hallmark)
  } else { empty_gsea_df() }
  for (ont in names(go_res)) { all_results[[paste0("GO_", ont)]] <- go_res[[ont]] }
  all_results[["KEGG"]] <- kegg_res
  all_results[["Reactome"]] <- reactome_res
  all_results[["Hallmark"]] <- hallmark_res
  all_results[["Column_Description"]] = description_gsea_df()
  write_xlsx(all_results, paste0(path, "3.GSEA/0.excel/", prefix, ".xlsx"))
  GSEA_result = all_results[[prefix_suffix]]
  draw_gsea_summary_plot(GSEA_result, prefix_suffix, path, Comparison_name, top_n = 30)
}

draw_Specific_GSEA_plot = function(path, data, Name, P.adjust = 0.05) {
  Genesets <- data$ID
  for (Geneset in Genesets) {
    try({
      Geneset_description = data$Description[data$ID == Geneset]
      NES = format(data$NES[data$ID == Geneset], digits = 4)
      FDR = format(data$p.adjust[data$ID == Geneset], scientific = TRUE, digits = 4)
      pvalue = format(data$pvalue[data$ID == Geneset], scientific = TRUE, digits = 4)
      setSize = data$setSize[data$ID == Geneset]
      color = ifelse(data$p.adjust[data$ID == Geneset] < P.adjust, "red", "black")
      make_folder(path, "3.GSEA/"); make_folder(path, "3.GSEA/2.GeneSet/")
      make_folder(path, paste0("3.GSEA/2.GeneSet/", Geneset_description))
      pdf(file = paste0(path, "3.GSEA/2.GeneSet/", Geneset_description, "/", Name, ".pdf"), width = 5, height = 5)
      p = gseaplot2(data, geneSetID = Geneset, subplots = 1:3, title = "", ES = "line")
      p[[1]] = p[[1]] +
        annotate("text", x = Inf, y = Inf, label = paste0("p-value: ", pvalue), hjust = 1, vjust = 2, size = 3) +
        annotate("text", x = Inf, y = Inf, label = paste0("FDR: ", FDR), hjust = 1, vjust = 4, size = 3, color = color) +
        annotate("text", x = Inf, y = Inf, label = paste0("NES: ", NES), hjust = 1, vjust = 6, size = 3) +
        labs(title = Name, subtitle = paste0(Geneset_description, "\n (setsize: ", setSize, ") "), y = "Running Enrichment Score") +
        theme(plot.title = element_text(hjust = 0.5, size = 14), plot.subtitle = element_text(hjust = 0.5, size = 12, color = color))
      print(p); dev.off()
    })
  }
}