# SYMBOL 이름을 가진 ranked vector를 ENTREZID 기준 vector로 변환하는 함수
convert_vector_to_entrez <- function(genes, orgdb, keytype) {
  original_ids <- names(genes)

  conversion <- bitr(
    original_ids,
    fromType = keytype,
    toType = "ENTREZID",
    OrgDb = orgdb
  )

  matched <- conversion[conversion[[keytype]] %in% original_ids, ]
  matched$log2FC <- genes[matched[[keytype]]]
  entrez_vector <- setNames(matched$log2FC, matched$ENTREZID)
  entrez_vector <- sort(entrez_vector, decreasing = TRUE)

  return(entrez_vector)
}

# GSEA 결과가 없을 때 반환할 빈 데이터프레임을 생성하는 함수
empty_gsea_df <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    setSize = character(),
    enrichmentScore = character(),
    NES = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    rank = character(),
    leading_edge = character(),
    core_enrichment = character(),
    stringsAsFactors = FALSE
  )
}

# GSEA 결과 각 컬럼의 의미를 설명하는 데이터프레임을 생성하는 함수
description_gsea_df <- function() {
    data.frame(
    Column = c(
        "Reference",
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust",
        "qvalue",
        "rank",
        "leading_edge",
        "core_enrichment"
    ),
    Meaning = c(
        "https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html",
        "Term ID",
        "Term name",
        "Number of genes in the term",
        "Raw enrichment score (ES); reflects how concentrated gene ranks are at the top/bottom of the ranked list",
        "Normalized Enrichment Score (ES normalized by set size and permutations)",
        "The p-value of the ES is calculated using a permutation test",
        "Adjusted p-value (Benjamini-Hochberg method)",
        "Adjusted p-value (Estimated false discovery rate (FDR))",
        "Rank in the ordered gene list where the enrichment score peak occurs",
        paste0(
        "Summary of the leading edge subset contributing most to the enrichment score:\n",
        " - tags: % of gene hits before (for +ES) or after (for -ES) the peak\n",
        " - list: % of ranked genes before (for +ES) or after (for -ES) the peak\n",
        " - signal: combined enrichment strength"
        ),
        "The actual genes contributing to the enrichment signal (leading edge core genes)"
        )
    )
}

# GSEA 결과를 p값과 gene set size 기준으로 필터링하는 함수
filter_gsea_result <- function(df, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
  result <- df %>%
    arrange(p.adjust) %>%
    filter(p.adjust < pvalue_cutoff, setSize >= min_set_size, setSize <= max_set_size)
  return(result)
}

# 항목 수와 라벨 길이에 따라 플롯 크기를 동적으로 조정하는 함수
adjust_plot_size <- function(top_df) {
  width <- 8
  height <- 8
  
  max_label_length <- max(nchar(as.character(top_df$Description)))
  counts = nrow(top_df)

  adjusted_width <- width + max_label_length * 0.1
  adjusted_height <- height - (20 - counts) * 0.3
  
  return(list(width = adjusted_width, height = adjusted_height))
}

# GSEA 결과에서 상위 pathway를 선택하고 시각화용 순서로 정렬하는 함수
select_top_GSEA_abs <- function(GSEA_result, top_n = 30) {
    Len = nrow(GSEA_result)
    GSEA_result %>%
        data.frame() %>%
        arrange(p.adjust) %>%
        head(n = min(top_n,Len)) %>%
        arrange(NES) -> top_pathway_data
    top_pathway_data$description = paste0(top_pathway_data$Description," (setsize: ",top_pathway_data$setSize,") ")
    factor(top_pathway_data$description, levels = top_pathway_data$description) -> top_pathway_data$description
    return(top_pathway_data)
}

# 요약용 GSEA bar plot을 PDF로 저장하는 함수
draw_gsea_summary_plot <- function(GSEA_result, prefix_suffix, path, Comparison_name, top_n = 30, pvalue_cutoff = 0.05, min_set_size = 10, max_set_size = 500) {
        set.seed(116)
        sort_df <- GSEA_result %>%
                arrange(p.adjust)
        sort_df = filter_gsea_result(sort_df, pvalue_cutoff = pvalue_cutoff, min_set_size = min_set_size, max_set_size = max_set_size)
        top_df = select_top_GSEA_abs(sort_df, top_n = top_n)

        plot_size <- adjust_plot_size(top_df)
        pdf(file= paste0(path,"3.GSEA/1.figure/",Comparison_name,"_",prefix_suffix,".pdf"), width=plot_size$width, height=plot_size$height)
        
        p = ggplot(top_df ) +
        geom_bar(aes(x = description, y = NES, fill = p.adjust),
                stat='identity') +
        scale_fill_continuous(low = "red", high = "blue") +
        coord_flip(ylim = c(-3.2, 3.2)) +
        scale_y_continuous(breaks = seq(-3,3,1)) +
        theme_bw(14) +
        labs(title=paste0(Comparison_name,"_",prefix_suffix),
                x = "GO Term") +  
        theme(plot.title = element_text(hjust = 0.5),
                axis.title.y = element_blank())

        print(p)
        dev.off()
}

### GSEA 함수

# GO BP/CC/MF에 대해 preranked GSEA를 수행하는 함수
perform_go_gsea <- function(genes, keytype, orgdb, Comparison_name, path) {
  ontologies <- c("BP", "CC", "MF")
  go_results <- list()

  for (ont in ontologies) {
    ego <- gseGO(
      gene          = genes,
      OrgDb         = orgdb,
      ont           = ont,
      keyType       = keytype,
      pAdjustMethod = "BH"
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      result <- as.data.frame(ego)
      result_sorted <- result[order(result$p.adjust), ]
    } else {
      result_sorted <- empty_gsea_df()
    }
    go_results[[ont]] <- result_sorted
  }

  return(go_results)
}

# ranked gene vector를 기반으로 KEGG GSEA를 수행하는 함수
perform_kegg_gsea <- function(genes, keytype, orgdb, Org_code) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)

  ekegg <- gseKEGG(
    gene = conv,
    organism = Org_code,
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH"
  )

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result <- result[order(result$p.adjust), ]
  } else {
    result <- empty_gsea_df()
  }
  return(result)
}

# ranked gene vector를 기반으로 Reactome GSEA를 수행하는 함수
perform_reactome_gsea <- function(genes, keytype, orgdb, organism) {
  conv <- convert_vector_to_entrez(genes, orgdb, keytype)

  erct <- gsePathway(
    gene          = conv,
    organism      = organism,
    pAdjustMethod = "BH"
  )
  if (!is.null(erct) && nrow(as.data.frame(erct)) > 0) {
    ekegg_readable <- setReadable(erct, OrgDb = orgdb, keyType = "ENTREZID")
    result <- as.data.frame(ekegg_readable)
    result_sorted <- result[order(result$p.adjust), ]
  } else {
    result_sorted <- empty_gsea_df()
  }

  return(result_sorted)
}

# 사용자 지정 GMT gene set에 대해 preranked GSEA를 수행하는 함수
perform_custom_geneset_gsea <- function(genes, Custom_geneset_gmt) {

  if (length(genes) == 0) return(empty_gsea_df())

  Custom_df <- read.gmt(Custom_geneset_gmt)

  gsea_res <- GSEA(
    geneList = genes,
    TERM2GENE = Custom_df,
    pAdjustMethod = "BH"
  )

  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    gsea_df <- as.data.frame(gsea_res)
    return(gsea_df[order(gsea_df$p.adjust), ])
  } else {
    return(empty_gsea_df())
  }
}

# 특정 gene set GMT에 대해 GSEA를 수행하고 개별 plot을 저장하는 함수
run_specific_gsea_analysis <- function(path, data, Comparison_name, specific_geneset, pvalue_cutoff = 1) {
    try({
        if (is.null(specific_geneset)) {
            return(invisible(NULL))
        }

        Genes <- setNames(data$log2FoldChange, data$gene_id)
        Genes <- sort(Genes, decreasing = TRUE)

        hallmark_df <- read.gmt(specific_geneset)

        gsea_res <- GSEA(
            geneList = Genes,
            TERM2GENE = hallmark_df,
            pAdjustMethod = "BH",
            pvalueCutoff = pvalue_cutoff
        )

        draw_Specific_GSEA_plot(path = path, data = gsea_res, Name = Comparison_name)
    })
}

# GSEA 전체 분석을 수행하고 결과 엑셀 및 대표 figure를 저장하는 함수
# [MODIFIED] Hallmark_path 인자 추가 — 하드코딩된 경로 대신 외부에서 전달받음
GSEA_analysis = function(path, data, Organism, Comparison_name, Keytype, Top_n = 20, Custom_geneset_gmt = NULL, prefix_suffix = "GO_BP", Hallmark_path = NULL) {
    make_folder(path,"3.GSEA/")
    make_folder(path,"3.GSEA/0.excel/")
    make_folder(path,"3.GSEA/1.figure/")

    if (Organism == "human"){
        OrgDb = org.Hs.eg.db
        Org_code = "hsa"
    } else if (Organism == "mouse") {
        OrgDb = org.Mm.eg.db
        Org_code = "mmu"
    }

    # Hallmark GMT: 외부 인자 우선, 없으면 스킵
    Hallmark = Hallmark_path
    
    all_results = list()

    Genes <- setNames(data$log2FoldChange, data$gene_id)
    Genes <- sort(Genes, decreasing = TRUE)
    prefix <- Comparison_name

    go_res <- perform_go_gsea(Genes, Keytype, OrgDb)
    kegg_res <- perform_kegg_gsea(Genes, Keytype, OrgDb, Org_code)
    reactome_res <- perform_reactome_gsea(Genes, Keytype, OrgDb, Organism)

    # Hallmark GSEA (경로가 있고 파일이 존재할 때만 수행)
    if (!is.null(Hallmark) && file.exists(Hallmark)) {
        hallmark_res <- perform_custom_geneset_gsea(Genes, Hallmark)
    } else {
        hallmark_res <- empty_gsea_df()
    }

    for (ont in names(go_res)) {
        GO_name = paste0("GO_", ont)
        all_results[[GO_name]] <- go_res[[ont]]
    }
    all_results[["KEGG"]] <- kegg_res
    all_results[["Reactome"]] <- reactome_res
    all_results[["Hallmark"]] <- hallmark_res
    all_results[["Column_Description"]] = description_gsea_df()
    write_xlsx(all_results, paste0(path,"3.GSEA/0.excel/",prefix,".xlsx"))

    GSEA_result = all_results[[prefix_suffix]]
    draw_gsea_summary_plot(
        GSEA_result = GSEA_result,
        prefix_suffix = prefix_suffix,
        path = path,
        Comparison_name = Comparison_name,
        top_n = 30
    )
}

# 특정 gene set 하나씩 GSEA enrichment plot을 저장하는 함수
draw_Specific_GSEA_plot = function(path, data, Name, P.adjust = 0.05){
    Genesets <- data$ID
    for (Geneset in Genesets){
        try({
            Geneset_description = data$Description[data$ID == Geneset]
            NES = format(data$NES[data$ID == Geneset], digits = 4)
            FDR = format(data$p.adjust[data$ID == Geneset], scientific = TRUE, digits = 4)
            pvalue = format(data$pvalue[data$ID == Geneset], scientific = TRUE, digits = 4)
            setSize = data$setSize[data$ID == Geneset]
            
            if (data$p.adjust[data$ID == Geneset] < P.adjust){
                color = "red"
            } else {
                color = "black"
            }
            
            make_folder(path, paste0("3.GSEA/"))
            make_folder(path, paste0("3.GSEA/2.GeneSet/"))
            make_folder(path, paste0("3.GSEA/2.GeneSet/", Geneset_description))
            pdf(file = paste0(path, "3.GSEA/2.GeneSet/", Geneset_description, "/", Name, ".pdf"), width = 5, height = 5)
            
            p = gseaplot2(data, geneSetID = Geneset, subplots = 1:3, title = "", ES = "line")
            
            p[[1]] = p[[1]] +
                annotate("text", x = Inf, y = Inf, label = paste0("p-value: ", pvalue), hjust = 1, vjust = 2, size = 3) + 
                annotate("text", x = Inf, y = Inf, label = paste0("FDR: ", FDR), hjust = 1, vjust = 4, size = 3, color = color) +
                annotate("text", x = Inf, y = Inf, label = paste0("NES: ", NES), hjust = 1, vjust = 6, size = 3) +
                labs(title = Name, subtitle = paste0(Geneset_description, "\n (setsize: ", setSize, ") "), y = "Running Enrichment Score") + 
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      plot.subtitle = element_text(hjust = 0.5, size = 12, color = color))
            
            print(p)
            dev.off()
        })
    }
}