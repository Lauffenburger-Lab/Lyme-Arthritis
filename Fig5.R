colnames(data_corrected) <- gsub("ApoB 100", "ApoB", colnames(data_corrected))

functional_groups <- c("IgG1", "IgG2","IgG3","IgA1", "ADCD")

for (fxn in functional_groups){
  #Plot 3: % of ADCD zeros in each fluid/phenotype combo
  q <- data_corrected[,which(grepl(fxn,colnames(data_corrected)))]
  antigens <- unlist(str_split(
    colnames(q), " "))[seq(2,length(unlist(
      str_split(colnames(q), " "))),2)]
  
  plot_zero_df <- data.frame(JF_refr = rep(NA, length(colnames(q))),
                             JF_resp = rep(NA, length(colnames(q))),
                             serum_refr = rep(NA, length(colnames(q))),
                             serum_resp = rep(NA, length(colnames(q))))
  rownames(plot_zero_df) <- colnames(q)
  
  p = 0
  for (sel_fluid in unique(fluid)){
    for (sel_phenotype in unique(phenotype)){
      p = p+1
      i = 0
      for (antigen in unique(antigens)){
        i = i + 1
        data_tmp <- q[which(phenotype == sel_phenotype & fluid == sel_fluid),i]
        count <- length(which(data_tmp < 0.1))/length(data_tmp)
        plot_zero_df[i,p] <- 1-count #to change to fraction WITH activity
      }
    }
  }
  
  plot_zero_df$antigen <- antigens
  plot_zero_df <- melt(plot_zero_df)
  plot_zero_df$antigen <- as.factor(plot_zero_df$antigen)
  plot_zero_df$fluid <- unlist(str_split(
    plot_zero_df$variable, "_"))[seq(1,length(unlist(
      str_split(plot_zero_df$variable, "_"))),2)]
  plot_zero_df$phenotype <- unlist(str_split(
    plot_zero_df$variable, "_"))[seq(2,length(unlist(
      str_split(plot_zero_df$variable, "_"))),2)]
  
  p3 <- plot_zero_df %>%
    mutate(antigen = fct_reorder(antigen, -value, .fun='median')) %>%
    ggplot(aes(x = antigen, y = value, 
               color = phenotype, group = variable, alpha = fluid)) + 
    geom_line(size = 1.5, aes(linetype = fluid)) +
    ylim(0,1) +
    scale_color_manual(values = c('#363795', '#41bb93')) +
    scale_alpha_discrete(range = c(0.7, 0.9)) +
    ylab(paste0("Fraction with ", fxn, " activity")) +
    scale_linetype_manual(name = "Linetype", 
                          values = c("JF" = 'solid', "serum" = 'dotted')) +
    xlab("Antigen") +
    theme_classic() +
    labs(title = fxn) +
    theme(axis.title.y = element_text(size = 20),
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 315, vjust = 0.1),
          plot.title = element_text(hjust=0.5, size = 20))
  
  print(p3)
  
  #Plot delta between JF & serum for each phenotype
  delta_df <- plot_zero_df
  delta_df_a <- delta_df[which(delta_df$variable %in% c("JF_refr", "JF_resp")),]
  delta_df_b <- delta_df[which(delta_df$variable %in% c("serum_refr", "serum_resp")),]
  delta_df <- data.frame(delta = delta_df_b$value/delta_df_a$value,
                         antigen1 = delta_df_a$antigen,
                         #to confirm correct subtraction, should be same as antigen 1
                         antigen2 = delta_df_b$antigen, 
                         phenotype = delta_df_a$phenotype)
  
  p3 <- delta_df %>%
    mutate(antigen1 = fct_reorder(antigen1, -delta, .fun='median')) %>%
    ggplot(aes(x = antigen1, y = delta, 
               color = phenotype, group = phenotype)) + 
    geom_point(size = 2) +
    scale_color_manual(values = c("#363795", "#41bb93")) +
    ylab(paste0("Ratio of serum to JF with ", fxn, " activity")) +
    xlab("Antigen") +
    theme_classic() +
    labs(title = fxn) +
    theme(axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 315, vjust = 0.1),
          plot.title = element_text(hjust=0.5, size = 10))
  
  print(p3)
  
}