setwd("")
df <- read.csv("fungi.csv")

df <- df[df$tool == "kraken2",] 

head(df)

ggplot(df, aes(as.character(sample), genus, fill = counts)) +
  geom_tile() +
  scale_fill_continuous(values = viridis(256), na.value = "white")

df$genus <- df$genus %>% str_trim("both")

dcb <- levels(factor(df$sample))
dfres <- data.frame("genus" = NA)
for (i in dcb) {
  dfpr <- df[df$sample == i,-c(1,4)]
  colnames(dfpr)[2] <- i
  
  dfres <- full_join(dfres, dfpr, by = "genus")
}

dfres <- dfres[-1,] 
rownames(dfres) <- dfres$genus
dfres <- as.matrix(dfres[,-1])
       
dfres[is.na(dfres)] <- 0

colpal1 = c('#028590', '#FF5733', '#FFC300', '#900C3F', '#581845', '#2E86AB', '#F4D35E', '#16697A')
colpal2 = gpar(start = '#028590', end = '#0c3b45', n = 256)

heatmaply(as.matrix(dfres), ColSideColors = colpal1, scale = "row") 

res.pca2 <- prcomp(t(dfres), scale = F)
fviz_pca_ind(res.pca2, repel = T)
