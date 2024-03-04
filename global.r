list.of.packages <- c(
  "stringr", "tidyr", "dplyr", "readr", "enrichR", "shiny", "shinydashboard",
  "ComplexHeatmap", "ggplot2", "DT"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages) > 0) {
  for (package in new.packages) {
    if (package == "ComplexHeatmap") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools")
      }
      devtools::install_github("jokergoo/ComplexHeatmap")
    } else {
      install.packages(package)
    }
  }
}



options(timeout = max(300, getOption("timeout")))

library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(enrichR)

tempo_de_timeout <- 300

file1 <- "https://www.dropbox.com/scl/fi/oikgx40lhatme27vtu311/lineplot.rds?rlkey=9m57na6a7wz3865gcu55mdz9k&dl=1"
destfile1 <- "lineplot.rds"

file2 <- "https://www.dropbox.com/scl/fi/cfc6mqxqdlendd4z0t58a/pathways_all.rds?rlkey=jqvc2w4qbnca23x5c3szuqz7g&dl=1"
destfile2 <- "pathways_all.rds"

file3 <- "https://www.dropbox.com/scl/fi/r78egopzy8vngwyya990k/Model.csv?rlkey=les5na781uq54dt0a9q7fd26b&dl=1"
destfile3 <- "Model.csv"

file4 <- "https://www.dropbox.com/scl/fi/pt0nykdnai9xeng6oc0ag/CIFS.txt?rlkey=i1yx4k18azz2mlrz8yk3qgplt&dl=1"
destfile4 <- "CIFS.txt"

baixar_arquivo <- function(url, destfile, timeout) {
  tryCatch(
    {
      download.file(url, destfile, mode = "wb", timeout = timeout)
      return(TRUE)
    },
    error = function(e) {
      return(FALSE)
    }
  )
}

if (baixar_arquivo(file1, destfile1, timeout = 300)) {
  lineplot <- readRDS(destfile1)
} else {
  stop("Falha ao baixar o arquivo")
}


# Baixar o arquivo
if (baixar_arquivo(file2, destfile2, timeout = 300)) {
  vias <- readRDS(destfile2)
} else {
  stop("Falha ao baixar o arquivo")
}

if (baixar_arquivo(file3, destfile3, timeout = 300)) {
  cell_annotation <- read_csv(destfile3)
} else {
  stop("Falha ao baixar o arquivo")
}

if (baixar_arquivo(file4, destfile4, timeout = 300)) {
  CIFS <- read.delim(destfile4, header = FALSE, col.names = "GENES")
} else {
  stop("Falha ao baixar o arquivo")
}

cell_annotation <- cell_annotation %>%
  select(OncotreeCode, CCLEName)

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
vias_nomes <- dbs[, 3]
dbs <- as.list(dbs[, 3])
names(dbs) <- vias_nomes

rm(vias_nomes)

Genes <- list(
  "All" = "All",
  "CIFS" = CIFS
)

Morpheus <- function(y) {
  malignos <- y %>% select(ends_with("Malignant"))
  linhagens2 <- y %>% select(-ends_with("Malignant"))

  if (nrow(linhagens2) == nrow(malignos)) {
    correlation_2 <- as.data.frame(cor(linhagens2, malignos, method = "spearman"))
    correlation_2 <- cbind("CCLEName" = rownames(correlation_2), correlation_2)
    correlation_2$CCLEName <- rownames(correlation_2)
    cell_annotation <- cell_annotation[which(cell_annotation$CCLEName %in% rownames(correlation_2)), ]
    correlation_2 <- merge(correlation_2, cell_annotation, by = "CCLEName", all.x = TRUE)
    rownames(correlation_2) <- correlation_2$CCLEName
    correlation_2 <- correlation_2 %>% select(-"CCLEName")
    correlation_2 <- correlation_2[, c("OncotreeCode", setdiff(names(correlation_2), "OncotreeCode"))]
    colnames(correlation_2)[1] <- "Subtype"
    gc()
    return(correlation_2)
  } else {
    print("linhagens2 e malignos têm comprimentos diferentes. A análise não pode ser realizada.")
    return(NULL)
  }
}

cor_p_value <- function(dataframe, selected_cancer = NULL) {
  resultados <- data.frame(
    Human_Cancer = character(),
    Cell_lines = character(),
    r_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )

  # Filtrando linhas com base nas escolhas do usuário
  if (!is.null(selected_cancer)) {
    malignos <- dataframe %>% select(ends_with("Malignant"), all_of(selected_cancer))
  } else {
    malignos <- dataframe %>% select(ends_with("Malignant"))
  }

  linhagens2 <- dataframe %>% select(-ends_with("Malignant"))
  malignos <- mutate_all(malignos, function(x) as.numeric(as.character(x)))
  linhagens2 <- mutate_all(linhagens2, function(x) as.numeric(as.character(x)))

  for (col_malignos in colnames(malignos)) {
    for (col_linhagens2 in colnames(linhagens2)) {
      cor_test <- cor.test(malignos[[col_malignos]], linhagens2[[col_linhagens2]])

      resultados <- rbind(
        resultados,
        data.frame(
          Human_Cancer = col_malignos,
          Cell_lines = col_linhagens2,
          r_value = cor_test$estimate,
          p_value = cor_test$p.value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  resultados <- resultados %>% arrange(desc(r_value))
  gc()
  return(resultados)
}

pathways <- function(lista, user) {
  nomes <- rownames(lista)
  aa <- c(user)
  analise <- as.data.frame(enrichr(nomes, aa))
  via <- user
  colnames(analise) <- str_replace_all(colnames(analise), paste(user, ".", sep = ""), "")
  grafico <- analise %>% filter(Adjusted.P.value < 0.05)
  grafico$Combined.Score <- as.numeric(as.character(grafico$Combined.Score))
  grafico$Term <- str_wrap(grafico$Term, width = 45)
  grafico <- grafico %>% arrange(desc(Combined.Score))
  grafico$Term <- factor(grafico$Term, levels = grafico$Term[order(grafico$Combined.Score)])
  gc()
  return(grafico)
}

dataframe <- function(data) {
  df <- data.frame(
    "Genes" = rownames(data),
    "log2FC" = data
  )
  return(df)
}
