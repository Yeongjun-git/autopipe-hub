# RNA-Seq 데이터 처리 관련 라이브러리
library(edgeR)       # edgeR DEG 분석
library(DESeq2)      # 차등 발현 분석
library(tximport)    # Salmon, Kallisto 등에서의 데이터 가져오기
library(RUVSeq)      # RNA-Seq 데이터 정규화
library(sva)         # 배치 효과 제거

# 생물학적 데이터 분석 관련 라이브러리
library(clusterProfiler) # 기능적 풍부도 분석
library(topGO)           # GO 분석
library(fgsea)           # GSEA 분석
library(ReactomePA)      # Reactome Pathway 분석
library(org.Mm.eg.db)    # 마우스 유전자 데이터베이스
library(org.Hs.eg.db)    # 인간 유전자 데이터베이스
library(biomaRt)         # Ensembl 데이터베이스 접근

# 시각화 관련 라이브러리
library(ggplot2)         # 데이터 시각화
library(ggrepel)         # ggplot에서 텍스트 레이블 정리
library(pheatmap)        # 히트맵 생성
library(enrichplot)      # 풍부도 분석 시각화
library(ggridges)        # 리지 플롯 생성
library(RColorBrewer)    # 색상 팔레트 제공

# 데이터 처리 및 분석 관련 라이브러리
library(dplyr)       # 데이터 조작 및 처리
library(tibble)      # 데이터 프레임 확장
library(readr)       # 빠른 데이터 읽기
library(readxl)      # 엑셀 파일 읽기
library(writexl)     # 엑셀 파일 쓰기
library(openxlsx)    # 엑셀 파일 쓰기

# 문자열 처리 관련 라이브러리
library(stringr)         # 문자열 처리

# CLI argument parsing
library(jsonlite)        # JSON 파싱 (comb_set 등)