### DEG enrichment test

# Hallmark 경로 resolve: HALLMARK_OVERRIDE > default
resolve_hallmark_path <- function(organism) {
  if (exists("HALLMARK_OVERRIDE", envir = .GlobalEnv)) {
    override <- get("HALLMARK_OVERRIDE", envir = .GlobalEnv)
    if (!is.null(override) && override != "" && file.exists(override)) {
      return(override)
    }
  }
  if (organism == "human") {
    return("/input/hallmark/h.all.v2025.1.Hs.symbols.gmt")
  } else if (organism == "mouse") {
    return("/input/hallmark/mh.all.v2025.1.Mm.symbols.gmt")
  }
  return("")
}

convert_symbol_to_entrez <- function(symbols, orgdb, keytype) {
  conversion <- bitr(symbols, fromType = keytype, toType = "ENTREZID", OrgDb = orgdb)
  return(conversion)
}

empty_enrichment_df <- function() {
  data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    stringsAsFactors = FALSE
  )
}

description_enrichment_df <- function() {
  data.frame(
    Column = c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"),
    Meaning = c(
      "Term ID", "Term name",
      "Ratio of DEGs annotated to this term (Intersected DEGs / QuerySize)",
      "Ratio of background genes annotated to this term (TermSize / Background gene set size)",
      "Raw enrichment p-value with Hypergeometric test",
      "Adjusted p-value (Benjamini-Hochberg method)",
      "Adjusted p-value (Estimated false discovery rate (FDR))",
      "List of DEGs in this term",
      "Number of matched DEGs with term"
    )
  )
}

perform_go_enrichment <- function(genes, keytype, orgdb) {
  ontologies <- c("BP", "CC", "MF")
  go_results <- list()
  for (ont in ontologies) {
    ego <- enrichGO(gene = genes, OrgDb = orgdb, ont = ont, keyType = keytype, pAdjustMethod = "BH", readable = TRUE)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      result <- as.data.frame(ego)
      result_sorted <- result[order(result$p.adjust), ]
    } else {
      result_sorted <- empty_enrichment_df()
    }
    go_results[[ont]] <- result_sorted
  }
  return(go_results)
}

perform_kegg_enrichment <- function(genes, keytype, orgdb, Org_code) {
  conv <- convert_symbol_to_entrez(genes, orgdb, keytype)
  if (nrow(conv) == 0) return(empty_enrichment_df())
  ekegg <- enrichKEGG(gene = unique(conv$ENTREZID), organism = Org_code, keyType = "ncbi-geneid", pAdjustMethod = "BH")
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result <- result[order(result$p.adjust), ]
  } else {
    result <- empty_enrichment_df()
  }
  return(result)
}

perform_reactome_enrichment <- function(genes, keytype, orgdb, organism) {
  conv <- convert_symbol_to_entrez(genes, orgdb, keytype)
  if (nrow(conv) == 0) return(empty_enrichment_df())
  erct <- enrichPathway(gene = unique(conv$ENTREZID), organism = organism, pAdjustMethod = "BH", readable = TRUE)
  if (!is.null(erct) && nrow(as.data.frame(erct)) > 0) {
    result <- as.data.frame(erct)
    result_sorted <- result[order(result$p.adjust), ]
  } else {
    result_sorted <- empty_enrichment_df()
  }
  return(result_sorted)
}

perform_custom_geneset_enrichment <- function(genes, keytype, Custom_geneset_gmt) {
  if (length(genes) == 0) return(empty_enrichment_df())
  Custom_geneset <- read.gmt(Custom_geneset_gmt)
  enrich_res <- enricher(gene = genes, TERM2GENE = Custom_geneset, pAdjustMethod = "BH")
  if (!is.null(enrich_res) && nrow(as.data.frame(enrich_res)) > 0) {
    enrich_df <- as.data.frame(enrich_res)
    return(enrich_df[order(enrich_df$p.adjust), ])
  } else {
    return(empty_enrichment_df())
  }
}

plot_DEG_enrichment_analysis <- function(df, Top_n, path, prefix, direction, genes, prefix_suffix = "GO_BP") {
  try({
    df$Intersection <- as.numeric(str_split(df$GeneRatio, "/", simplify = TRUE)[, 1])
    df$QuerySize <- as.numeric(str_split(df$GeneRatio, "/", simplify = TRUE)[, 2])
    df$TermSize <- as.numeric(str_split(df$BgRatio, "/", simplify = TRUE)[, 1])
    df$generatio <- df$Intersection / df$QuerySize
    max <- min(Top_n, nrow(df))
    topn_df <- df[order(df$p.adjust), ][1:max, ]
    plot_size <- adjust_plot_size(topn_df)
    pdf(file = paste0(path, "2.DEG_enrichment/1.figure/", prefix, "_", prefix_suffix, ".pdf"),
        width = plot_size$width + 5, height = plot_size$height)
    p <- ggplot(topn_df, aes(x = generatio, 
                              y = reorder(paste0(Description, " (size: ", TermSize, ") "), -log10(p.adjust)), 
                              size = Count, color = -log10(p.adjust))) +
      geom_point(alpha = 0.6) + scale_color_gradient(low = "blue", high = "red") +
      theme_minimal() +
      labs(title = paste0(prefix, "_", prefix_suffix), x = "Intersection Size / Term Size", 
           y = "Term", size = "Intersection Size", color = "-Log p-value") +
      annotate("text", x = Inf, y = Inf, label = paste(direction, " regulated DEGs: ", length(genes)), vjust = 2, hjust = 1, size = 3) +
      theme(axis.text = element_text(size = 16), plot.title = element_text(hjust = 0.5))
    print(p)
    dev.off()
  })
}

GO_DEG_analysis = function(path, data, Organism, Comparison_name, Keytype, Top_n = 40, Custom_geneset_gmt = NULL) {
    make_folder(path, "2.DEG_enrichment/")
    make_folder(path, "2.DEG_enrichment/0.excel/")
    make_folder(path, "2.DEG_enrichment/1.figure/")
    if (Organism == "human") {
        OrgDb = org.Hs.eg.db; Org_code = "hsa"
    } else if (Organism == "mouse") {
        OrgDb = org.Mm.eg.db; Org_code = "mmu"
    }
    Hallmark <- resolve_hallmark_path(Organism)

    sets = c("UP", "DOWN")
    for (direction in sets) {
        all_results <- list()
        genes <- data$gene_id[data$diffexpressed == direction]
        prefix <- paste(Comparison_name, direction, sep = "__")
        go_res <- perform_go_enrichment(genes, Keytype, OrgDb)
        kegg_res <- perform_kegg_enrichment(genes, Keytype, OrgDb, Org_code)
        reactome_res <- perform_reactome_enrichment(genes, Keytype, OrgDb, Organism)
        hallmark_res <- if (Hallmark != "" && file.exists(Hallmark)) {
          perform_custom_geneset_enrichment(genes, Keytype, Hallmark)
        } else { empty_enrichment_df() }
        for (ont in names(go_res)) {
            all_results[[paste0("GO_", ont)]] <- go_res[[ont]]
        }
        all_results[["KEGG"]] <- kegg_res
        all_results[["Reactome"]] <- reactome_res
        all_results[["Hallmark"]] <- hallmark_res
        if (!is.null(Custom_geneset_gmt)) {
            all_results[["Custom"]] <- perform_custom_geneset_enrichment(genes, Keytype, Custom_geneset_gmt)
        }
        all_results[["Column_Description"]] = description_enrichment_df()
        write_xlsx(all_results, paste0(path, "2.DEG_enrichment/0.excel/", prefix, ".xlsx"))
        prefix_suffix = "GO_BP"
        plot_DEG_enrichment_analysis(all_results[[prefix_suffix]], Top_n, path, prefix, direction, genes, prefix_suffix)
    }
}