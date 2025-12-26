#!/opt/conda/bin/Rscript
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))  
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
library(optparse)
option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = FALSE,
                help = "input file path"),
    make_option(c("-o", "--output"), type = "character", default = FALSE,
              help = "output file path"),
    make_option(c("-p", "--prefix"), type = "character", default = FALSE,
              help = "prefix of output file"),
    make_option(c("-t", "--type"), type = "character", default = FALSE,
              help = "type of mod (m6A, m5C, AtoI, psi, etc.)")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input  <-  opt$input
output <- opt$output
prefix <- opt$prefix
type <- opt$type

# library packages
library(dplyr)
library(purrr)
library(scales)
library(ggplot2)
library(magrittr)
library(legendBaseModel)

plot_metagene_Rd <- function(input_data, output_path, prefix, mod_type = "m6A") {
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    tryCatch({
    # START
    cat(crayon::green(paste0(prefix, " Start ploting Metagene\n")))
    colour_fill <- c("black", "red")
    names(colour_fill) <- input_data$group  %>% unique()
    p1 <- ggplot(input_data, aes(x=value))+
        geom_line(aes(colour = group), stat = "density", adjust = 2) +
        scale_colour_manual(values=colour_fill) +
        geom_vline(xintercept = 1:2, col = "red", linetype="dashed") + 
        annotate("text", x = 0.5, y = -0.14, label = "5'UTR") +
        annotate("text", x = 1.5, y = -0.14, label = "CDS") +
        annotate("text", x = 2.5, y = -0.14, label = "3'UTR")  +
        annotate("rect", xmin = 0, xmax = 1, ymin = -0.08, ymax = -0.04, alpha = .99, colour = "black") +
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.08, ymax = -0.04, alpha = .99, colour = "black") +
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = 0, alpha = .2, colour = "black") +
        xlab(paste0(mod_type, " metagene")) +
        ylab("Frequency") +
        guides(colour = guide_legend(title = NULL)) +
        theme_bw() +
        theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
        theme(legend.background = element_blank()) 
    ggsave(paste0(output_path, "/", prefix, "_", mod_type, "_metagene.pdf"), p1, width = 8, height = 6)
    ggsave(paste0(output_path, "/", prefix, "_", mod_type, "_metagene.png"), p1, width = 8, height = 6)
    ggsave(paste0(output_path, "/", prefix, "_", mod_type, "_metagene.tiff"), p1, width = 8, height = 6)

    cat(crayon::green(paste0(prefix, " Metagene analysis done\n")))
    # END
    return("sucess")},
    error = function(e){print(e$message);message(return("failled"))})
}
restructure_coord <- function(input_path, group, type) {
    if (type == "m6A") {
        m6a.dist <- read.delim(input_path, header = T) %>%
            dplyr::filter(utr5_size > utr3_size) %>% 
            dplyr::mutate(
                trx_len = utr5_size + cds_size + utr3_size)
    } else if (type == "m5C") {
        m6a.dist <- read.delim(input_path, header = T) %>%
            dplyr::filter(cds_size > utr3_size & utr3_size > utr5_size) %>% 
            dplyr::mutate(
                trx_len = utr5_size + cds_size + utr3_size)
    } else if (type == "psi") {
        m6a.dist <- read.delim(input_path, header = T) %>%
            dplyr::filter(cds_size > utr3_size & utr3_size > utr5_size) %>% 
            dplyr::mutate(
                trx_len = utr5_size + cds_size + utr3_size)
    } else if (type == "AtoI") {
        m6a.dist <- read.delim(input_path, header = T) %>% 
            dplyr::filter(cds_size > utr3_size & utr3_size > utr5_size) %>% 
            dplyr::mutate(
                trx_len = utr5_size + cds_size + utr3_size)
    }

    temp.df <- m6a.dist %>%
        dplyr::select(c("gene_name", "refseqID", "trx_len")) %>%
        dplyr::rename(setNames(colnames(.), c("gene_name", "gid", "trx_len"))) %$%
        .[order(.$gene_name, .$gid, -.$trx_len), ] %$%
        .[!duplicated(.$gene_name), ]

    m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid, ]
    utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    # print(paste0("utr5.SF: ", utr5.SF, " utr3.SF: ", utr3.SF))
    utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
    cds.m6a.dist <- m6a.dist[m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
    utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2 & m6a.dist$rel_location < 3, ]
    utr5.m6a.dist$rel_location <- scales::rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0, 1))
    utr3.m6a.dist$rel_location <- scales::rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2, 3))
    m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)

    m6a_data <- data.frame("value"= m6a.metagene.coord) %>% 
        dplyr::mutate("group"= group)
    return(m6a_data)
}
auto_meta_plot <- function(input_path, output_path, prefix, mod_type = "m6A") {
    # test
    # 2. check output path
    options(warn = -1)
    ifelse(dir.exists(output_path), "Dir exist alreadly!", dir.create(output_path, recursive = TRUE))
    tryCatch({
    # START
    sample_all <- legendBaseModel::make_sample_sheet(input_path, paste0("_", mod_type, "_abs_dist.txt")) %>%
        dplyr::mutate(
            group = sample,
            type = mod_type,
            dist_measures = purrr::pmap(list(sample_path, group, type), restructure_coord)
        ) %>% 
        dplyr::mutate(
            output_plot = output_path,
            prefix_plot = paste0(prefix, "_", group),
            plot_meta <- purrr::pmap(list(dist_measures, output_plot, prefix_plot, type), plot_metagene_Rd)
        )
    # group plot
    legendBaseModel::rbind_merge(sample_all, "dist_measures") %>% 
        plot_metagene_Rd(., output_path, paste0(prefix, "_group"))
    # END
    return("sucess")},
    error = function(e){print(e$message);message(return("failled"))})
}


# main run
auto_meta_plot(input, output, prefix, type)