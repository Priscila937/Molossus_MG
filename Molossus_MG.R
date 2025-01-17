#carregar os pacotes
library(geomorph)
library(Morpho)
library(ggplot2)
library(phytools)
library(RColorBrewer)
install.packages("ape")
library(ape)
install.packages("dendextend")
library(dendextend)
install.packages("tidytree")
library(treeio)
install.packages("ggplotify")
library(ggtree)

#mudar o diretório de trabalho
setwd("")

#ler o arquivo de landmarks
Especies_Mol_I<-readland.tps("Especies_Molossus_inferior.tps", readcurves = F)
Especies_Mol_I.sobrepo<-procSym(Especies_Mol_I) #função do pacote Morpho
Especies_Mol_I.sobrepo2<-gpagen(Especies_Mol_I) #geomorph

plotAllSpecimens(Especies_Mol_I)
plot(Especies_Mol_I.sobrepo2)
dim(Especies_Mol_I)

sizeEI<-Especies_Mol_I.sobrepo$size #ver escala

#categorizar os grupos
gruposICVA<-as.factor(c(rep("Mcur", 2), rep("Mmil", 3), rep("Mpar", 5),rep("Mmel", 3), 
                        rep("M.flu", 6), rep("M.azt", 12), rep("M.coi", 3), rep("M.pre", 7), 
                        rep("M.mol", 174), rep("M.ruf", 84)))

PCA.gls <- gm.prcomp(Especies_Mol_I.sobrepo$rotated)
scoresIR<-Especies_Mol_I.sobrepo$PCscores[,1:12]

#Separando valores de PC1 e PC2 para todos os especimes
PC1e2<-Especies_Mol_I.sobrepo$PCscores[,1:2]

#Estabelecendo cores para as espécies
cores_especies <- c(brewer.pal(n = 9, name = "Set1"), "black") 
names(cores_especies) <- c("Mcur", "Mmil", "Mpar", "Mmel",
                           "M.flu", "M.azt", "M.coi", "M.pre", "M.mol", "M.ruf")

# Definir as 10 formas diferentes para as espécies
simbolos_especies <- c(
  Mcur = 0,  # Triângulo
  Mmil = 1,  # Círculo
  Mpar = 2,  # Quadrado
  Mmel = 3,  # Diamante
  M.flu = 4,  # Triângulo invertido
  M.azt = 5,  # Mais
  M.coi = 6,  # Cruz
  M.pre = 7,  # Estrela
  M.mol = 8, # Asterisco
  M.ruf = 9   # Triângulo à direita
)

#Plotar o grafico da PCA
# Gerar o gráfico
ggplot(PC1e2, aes(x = PC1, y = PC2, color = gruposICVA, shape = gruposICVA)) +
  geom_point(aes(color = gruposICVA, shape = gruposICVA), size = 5, stroke = 1.5, na.rm = TRUE) +
  scale_shape_manual(values = simbolos_especies) + 
  scale_color_manual(values = cores_especies) +
  labs(x = "PCA1", y = "PCA2") +
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Aumenta tamanho dos símbolos na legenda
    legend.text = element_text(size = 14),  # Texto maior na legenda
    axis.text = element_text(size = 12),  # Texto dos eixos maior
    axis.title = element_text(size = 16)  # Título dos eixos maior
  )


# Criar um data frame com os dados
dados <- data.frame(
  Grupos = gruposICVA,
  TamanhoCentroide = sizeEI
)

# Gerar gráfico box plot do tipo violin
violin_plot <- ggplot(dados, aes(x = Grupos, y = TamanhoCentroide, fill = Grupos)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = cores_especies) +
  theme_minimal() +
  labs(
    x = "cores_especies",
    y = "Tamanho do Centroide"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

print(violin_plot)

#CVA entre os grupos com teste de permutacao de 1000 vezes e cross-validation)            
cvall<-CVA(scoresIR,gruposICVA,rounds=1000,cv=T) 
typprobs <- typprobClass(cvall$CVscores,groups=gruposICVA)
print(typprobs)

# plot CVA
#CVA com ggplot
cva.scores <- as.data.frame(cvall$CVscores) #transformar a matriz em dataframe
#cvall é o resultado da análise de CVA de qualquer coisa
#cvall$CVscores será usado para extrair somente os eixos da CVA

cva.scores$species<-gruposICVA #adicionar os grupos no dataframe

# Criar o gráfico com os dados da CVA
ggplot(cva.scores, aes(x = `CV 1`, y = `CV 2`)) +
  geom_point(aes(color = species, shape = species), size = 5, stroke = 1.5, na.rm = TRUE) +
  scale_shape_manual(values = simbolos_especies) +
  scale_color_manual(values = cores_especies)+
  coord_fixed(ratio = 1) +
  labs(x = "CV 1", y = "CV 2", color = "Espécies", shape = "Espécies") +
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Aumenta tamanho dos símbolos na legenda
    legend.text = element_text(size = 14),  # Texto maior na legenda
    axis.text = element_text(size = 12),  # Texto dos eixos maior
    axis.title = element_text(size = 16)  # Título dos eixos maior
  )


#ANOVA testanto o dimorfismo
gdf.inferior<-geomorph.data.frame(coords=scoresIR, grupos=gruposICVA, size=sizeEI)

fit1.inferior <- procD.lm(coords ~ size, 
                          data = gdf.inferior, iter = 9999, 
                          RRPP = F, print.progress = F) # randomize raw values
summary(fit1.inferior)


#Testando a significancia do centroide
fit2.inferior <- procD.lm(size ~ grupos, 
                          data = gdf.inferior, iter = 9999, 
                          RRPP = F, print.progress = F)

summary(fit2.inferior)

#extraindo os resíduos
resI<-fit1.inferior$residuals
cvallCR<-CVA(resI, gruposICVA, rounds=1000, cv=T) #com residuos

#CVA com resíduos
cva.scoresCR <- as.data.frame(cvallCR$CVscores) #transformar a matriz em dataframe
cva.scoresCR$species<-gruposICVA #adicionar os grupos no dataframe

# Criar o gráfico com os dados da CVA com resíduos
ggplot(cva.scoresCR, aes(x = `CV 1`, y = `CV 2`)) + 
  geom_point(aes(color = species, shape = species), size = 5, stroke = 1.5) + 
  scale_shape_manual(values = simbolos_especies) + 
  scale_color_manual(values = cores_especies)+
  coord_fixed(ratio = 1) + 
  labs(x = "CV 1", y = "CV 2", color = "Espécies", shape = "Espécies") + 
  theme_minimal() +
  theme(
    legend.key.size = unit(1.5, "cm"),  # Aumenta tamanho dos símbolos na legenda
    legend.text = element_text(size = 14),  # Texto maior na legenda
    axis.text = element_text(size = 12),  # Texto dos eixos maior
    axis.title = element_text(size = 16)  # Título dos eixos maior
  )
  
#analise de integração das vistas do cranio
two.b.pls(scoresIR, scoresSR, iter = 999, seed = NULL, print.progress = TRUE)

# Escores de regressão
M <- Especies_Mol_I.sobrepo2$consensus
plethAllometry <- procD.lm(coords ~ grupos, data=gdf.inferior)

# Criar o data frame de referência para espécies
especies_info <- data.frame(
  grupo = names(cores_especies),  # Nomes das espécies
  cor = unname(cores_especies),   # Cores associadas
  shape = unname(simbolos_especies)  # Shapes associados
)

# Mapear cores e shapes com base nos gruposICVA
cores_plot <- especies_info$cor[match(gruposICVA, especies_info$grupo)]
shapes_plot <- especies_info$shape[match(gruposICVA, especies_info$grupo)]

# Verificar se o mapeamento foi feito corretamente
head(data.frame(gruposICVA, cores_plot, shapes_plot))

allom.plot <- plot(
  plethAllometry,
  type = "regression",
  predictor = log(gdf.inferior$size),
  reg.type = "RegScore",
  pch = shapes_plot,  # Shapes
  col = cores_plot, 
  lwd = 2             # Espessura da linha
)

#analise da forma para os grupos de espécies. Espécies menores versus maiores  
PC <- PCA.gls$x[,1]
preds <- shape.predictor(Especies_Mol_I.sobrepo2$coords, x= cvall$CVscores[,1], Intercept = FALSE, 
                         pred1 = min(cvall$CVscores[,1]), pred2 = max(cvall$CVscores[,1])) # PC 1 extremes, more technically
plotRefToTarget(M, preds$pred1, mag=1)
plotRefToTarget(M, preds$pred2, mag=1)

#plotar pontos de referencia para o preds1
plotRefToTarget(M, preds$pred1, mag=1)
segments(preds$pred1[1, 1], preds$pred1[1, 2], preds$pred1[2, 1], preds$pred1[2, 2], lwd = 5, col = "blue")  # 1 to 2  # 1 to 2
segments(preds$pred1[2, 1], preds$pred1[2, 2], preds$pred1[7, 1], preds$pred1[7, 2], lwd = 5, col = "blue")  # 2 to 7
segments(preds$pred1[7, 1], preds$pred1[7, 2], preds$pred1[11, 1], preds$pred1[11, 2], lwd = 5, col = "blue")  # 7 to 11
segments(preds$pred1[11, 1], preds$pred1[11, 2], preds$pred1[17, 1], preds$pred1[17, 2], lwd = 5, col = "blue")  # 11 to 17
segments(preds$pred1[17, 1], preds$pred1[17, 2], preds$pred1[15, 1], preds$pred1[15, 2], lwd = 5, col = "blue")  # 17 to 15
segments(preds$pred1[15, 1], preds$pred1[15, 2], preds$pred1[16, 1], preds$pred1[16, 2], lwd = 5, col = "blue")  # 15 to 16
segments(preds$pred1[16, 1], preds$pred1[16, 2], preds$pred1[14, 1], preds$pred1[14, 2], lwd = 5, col = "blue")  # 16 to 14
segments(preds$pred1[14, 1], preds$pred1[14, 2], preds$pred1[18, 1], preds$pred1[18, 2], lwd = 5, col = "blue")  # 14 to 18
segments(preds$pred1[18, 1], preds$pred1[18, 2], preds$pred1[12, 1], preds$pred1[12, 2], lwd = 5, col = "blue")  # 18 to 12
segments(preds$pred1[12, 1], preds$pred1[12, 2], preds$pred1[8, 1], preds$pred1[8, 2], lwd = 5, col = "blue")  # 12 to 8
segments(preds$pred1[8, 1], preds$pred1[8, 2], preds$pred1[3, 1], preds$pred1[3, 2], lwd = 5, col = "blue")  # 8 to 3
segments(preds$pred1[3, 1], preds$pred1[3, 2], preds$pred1[1, 1], preds$pred1[1, 2], lwd = 5, col = "blue")  # 3 to 1
segments(preds$pred1[15, 1], preds$pred1[15, 2], preds$pred1[13, 1], preds$pred1[13, 2], lwd = 5, col = "blue")  # 15 to 13
segments(preds$pred1[13, 1], preds$pred1[13, 2], preds$pred1[14, 1], preds$pred1[14, 2], lwd = 5, col = "blue")  # 13 to 14
segments(preds$pred1[4, 1], preds$pred1[4, 2], preds$pred1[5, 1], preds$pred1[5, 2], lwd = 5, col = "blue")  # 4 to 5
segments(preds$pred1[5, 1], preds$pred1[5, 2], preds$pred1[6, 1], preds$pred1[6, 2], lwd = 5, col = "blue")  # 5 to 6
segments(preds$pred1[4, 1], preds$pred1[4, 2], preds$pred1[9, 1], preds$pred1[9, 2], lwd = 5, col = "blue")  # 4 to 9
segments(preds$pred1[6, 1], preds$pred1[6, 2], preds$pred1[10, 1], preds$pred1[10, 2], lwd = 5, col = "blue")  # 6 to 10

#plotar pontos de referencia para o preds2
plotRefToTarget(M, preds$pred2, mag=1)
segments(preds$pred2[1, 1], preds$pred2[1, 2], preds$pred2[2, 1], preds$pred2[2, 2], lwd = 5, col = "red")  # 1 to 2  # 1 to 2
segments(preds$pred2[2, 1], preds$pred2[2, 2], preds$pred2[7, 1], preds$pred2[7, 2], lwd = 5, col = "red")  # 2 to 7
segments(preds$pred2[7, 1], preds$pred2[7, 2], preds$pred2[11, 1], preds$pred2[11, 2], lwd = 5, col = "red")  # 7 to 11
segments(preds$pred2[11, 1], preds$pred2[11, 2], preds$pred2[17, 1], preds$pred2[17, 2], lwd = 5, col = "red")  # 11 to 17
segments(preds$pred2[17, 1], preds$pred2[17, 2], preds$pred2[15, 1], preds$pred2[15, 2], lwd = 5, col = "red")  # 17 to 15
segments(preds$pred2[15, 1], preds$pred2[15, 2], preds$pred2[16, 1], preds$pred2[16, 2], lwd = 5, col = "red")  # 15 to 16
segments(preds$pred2[16, 1], preds$pred2[16, 2], preds$pred2[14, 1], preds$pred2[14, 2], lwd = 5, col = "red")  # 16 to 14
segments(preds$pred2[14, 1], preds$pred2[14, 2], preds$pred2[18, 1], preds$pred2[18, 2], lwd = 5, col = "red")  # 14 to 18
segments(preds$pred2[18, 1], preds$pred2[18, 2], preds$pred2[12, 1], preds$pred2[12, 2], lwd = 5, col = "red")  # 18 to 12
segments(preds$pred2[12, 1], preds$pred2[12, 2], preds$pred2[8, 1], preds$pred2[8, 2], lwd = 5, col = "red")  # 12 to 8
segments(preds$pred2[8, 1], preds$pred2[8, 2], preds$pred2[3, 1], preds$pred2[3, 2], lwd = 5, col = "red")  # 8 to 3
segments(preds$pred2[3, 1], preds$pred2[3, 2], preds$pred2[1, 1], preds$pred2[1, 2], lwd = 5, col = "red")  # 3 to 1
segments(preds$pred2[15, 1], preds$pred2[15, 2], preds$pred2[13, 1], preds$pred2[13, 2], lwd = 5, col = "red")  # 15 to 13
segments(preds$pred2[13, 1], preds$pred2[13, 2], preds$pred2[14, 1], preds$pred2[14, 2], lwd = 5, col = "red")  # 13 to 14
segments(preds$pred2[4, 1], preds$pred2[4, 2], preds$pred2[5, 1], preds$pred2[5, 2], lwd = 5, col = "red")  # 4 to 5
segments(preds$pred2[5, 1], preds$pred2[5, 2], preds$pred2[6, 1], preds$pred2[6, 2], lwd = 5, col = "red")  # 5 to 6
segments(preds$pred2[4, 1], preds$pred2[4, 2], preds$pred2[9, 1], preds$pred2[9, 2], lwd = 5, col = "red")  # 4 to 9
segments(preds$pred2[6, 1], preds$pred2[6, 2], preds$pred2[10, 1], preds$pred2[10, 2], lwd = 5, col = "red")  # 6 to 10

#Verificar a distancia de Mahalanobis
# Plot distâncias de Mahalahobis Cvall
dendroS=hclust(cvall$Dist$GroupdistMaha)
dendroS$labels=levels(gruposICVA)
par(mar=c(4,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="Espécies",
     ylab='Mahalahobis distance')

# plot distâncias de Mahalahobis cvall resíduos
dendroS=hclust(cvallCR$Dist$GroupdistMaha)
dendroS$labels=levels(gruposICVA)
par(mar=c(4,4.5,1,1))
dendroS=as.dendrogram(dendroS)
plot(dendroS, main='',sub='', xlab="Espécies",
     ylab='Mahalahobis distance')

#install.packages("dendextend") para comparar os dendrogramas anteriores
#library(dendextend)
# Gerar o primeiro dendrograma, com efeito de tamanho
dendroS1 <- hclust(cvall$Dist$GroupdistMaha)
dendroS1$labels <- levels(gruposICVA)
dendroS1 <- as.dendrogram(dendroS1)

# Gerar o segundo dendrograma, sem efeito de tamanho
dendroS2 <- hclust(cvallCR$Dist$GroupdistMaha)
dendroS2$labels <- levels(gruposICVA)
dendroS2 <- as.dendrogram(dendroS2)

# Comparar os dois dendrogramas
dend_list <- dendlist(dendroS1, dendroS2)

# Plotar os dois dendrogramas lado a lado com linhas de ligação
tanglegram(dend_list, 
           highlight_distinct_edges = FALSE, 
           common_subtrees_color_lines = TRUE, 
           common_subtrees_color_branches = TRUE, 
           #main_left = "Dendrograma com efeito de tamanho",
           #main_right = "Dendrograma sem efeito de tamanho",
           lab.cex = 1.2, edge.lwd = 2)

#ler a árvore filogenética
# Instalar e carregar o pacote treeio e o ggplotify, se ainda não o tiver
# Carregar a árvore
#Adicionar o caminho no seu pc para chegar ao arquivo
tree <- read.beast("")

# Verificar a estrutura da árvore
str(tree)

# Conferir as idades dos nós
tree@phylo$node.label
tip <- c("MH185125_MH410726_Malvarezi", "EF080483_MH058051_Mfentoni", "MG191810_KM387368_Mbondae",
         "MH185179_MH058094_Msinaloae", "MH185186_MH058091_P.centralis", "KX355065_MH058057_Mverrilli", "JF454657_MH058046_E.auripendulus")

#excluir os táxons não amostrados no conjunto de dados morfológico
batstotal.tre <- drop.tip(tree, tip)

ggtree(batstotal.tre) +
  geom_tiplab(align = TRUE, linetype = 'dashed', linesize = .3)

#boxplot de uma característica alvo
df<-data.frame(x=gruposICVA, y=sizeEI)

#Calculando as médias dos centroides
sizecur<-mean(df[1:2,2])
sizemil<-mean(df[3:5,2])
sizepar<-mean(df[6:10,2])
sizemel<-mean(df[11:13,2])
sizeflu<-mean(df[14:19,2])
sizeazt<-mean(df[20:31,2])
sizecoi<-mean(df[32:34,2])
sizepre<-mean(df[35:41,2])
sizemol<-mean(df[42:215,2])
sizeruf<-mean(df[216:299,2])  

#Calcular as medias do centroide para cada especie e associar o valor a cada uma
mediascent<-c(sizemil, sizepar, sizeazt,  sizemel, sizecoi, sizeflu, sizeruf, sizecur, sizepre, sizemol)
especiesmedia<-c("MH185146_MH058056_Mmilleri_Cuba", 
                 "129_Mparanaensis", 
                 "MH185133_MH058047_Maztecus_Mexico", 
                 "BDN1_Mmelini_BrasilPR",
                 "CA28_MH058078_Mcoibensis",  
                 "AMA114_Mfluminensis", 
                 "CESC64_Mrufus", 
                 "MH185139_MH058050_Mcurrentium_Paraguay",  
                "MH185168_MH058079_Mpretiosus_Nicaragua", 
                "CUMA73_Mmolossus_BrazilMA")

names(mediascent)<-especiesmedia
# Converter a árvore 'treedata' para 'phylo'
bats.tree <- as.phylo(batstotal.tre)
#verificar se a arvore e da classe phylo
class(bats.tree)

#Calculando o estado ancestral
ancestral_states <- fastAnc(bats.tree, mediascent, vars=FALSE, CI=FALSE)

# Visualização simples com phytools do estado ancestral
plot(bats.tree)
nodelabels(ancestral_states, frame="n", cex=0.8)

## plot contMap
obj<-contMap(bats.tree, mediascent)
ln.size<-log(setNames(mediascent,
                          rownames(especiesmedia)))

names(ln.size)<-especiesmedia
mammal.contMap<-contMap(bats.tree,
                        ln.size,plot=FALSE,res=200)

## change color scheme
mammal.contMap<-setMap(mammal.contMap,
                       c("white","#FFFFB2","#FECC5C","#FD8D3C",
                         "#E31A1C"))
plot(mammal.contMap,fsize=c(0.7,0.8),
     leg.txt="log(centroid size)")
par(mar=c(5.1,4.1,4.1,2.1))

#Testando sinal filogenetico
ps<-physignal(mediascent, bats.tree,iter=9999)
ps

# Estimar o tamanho do efeito
ps_z <- physignal.z(mediascent, bats.tree)
ps_z

#Ler o arquivo com os dados de PC
PC<- read.table("Especies-PC1-PC2.txt", header = TRUE)

# Verificar se todos os nomes de táxons em PC estão presentes em tree$tip.label
if (!all(PC$species %in% bats.tree$tip.label)) {
  stop("Há inconsistências entre os nomes dos táxons nos dados de PC e na árvore filogenética.")
}

## Definir a coluna 1 como os nomes das linhas do data frame
rownames(PC) <- PC$Especies

# Reordenar os dados de PC para corresponder à ordem dos rótulos dos táxons na árvore
PC2 <- PC[match(bats.tree$tip.label, PC$Especies), ]

#phylomorphospace
phylomorphospace(bats.tree, PC[, 2:3], colors=setNames(c("blue","red"),c(0,1)))

#Comparing net rates of evolution among traits on phylogenies
#Criando grupo
gps <- c("129_Mparanaensis",
         "AMA114_Mfluminensis",
         "BDN1_Mmelini_BrasilPR",
         "CA28_MH058078_Mcoibensis",
         "CESC64_Mrufus",
         "CUMA73_Mmolossus_BrazilMA",
         "MH185133_MH058047_Maztecus_Mexico",
         "MH185139_MH058050_Mcurrentium_Paraguay",
         "MH185146_MH058056_Mmilleri_Cuba",
         "MH185168_MH058079_Mpretiosus_Nicaragua"
)

#verificar se a arvore e da classe phylo
class(bats.tree)
bats.tree <- as.phylo(bats.tree) 
bats.tree$tip.label #Validando os nomes das pontas da árvore
names(PC2) #confirmando o nome dos vetores


PCs <- as.matrix(PC2[, 2:3])  # Selecionar apenas PC1 e PC2 e transformar em matriz
str(PCs)# Verificar se PCs agora é uma matriz
rownames(PCs) <- PC2$Especies  # Atribuir os nomes das espécies como nomes das linhas
rownames(PCs) <- PC2$Especies  # Reatribua os nomes de linha, se necessário
rownames(PCs)  # Confirmar que os nomes das linhas correspondem às espécies
all(bats.tree$tip.label %in% rownames(PCs))  # Verificar se todos os tips da árvore estão no PCs
all(bats.tree$tip.label %in% gps) 
all(gps %in% PCs) #Comparar os nomes mediascent com gps se esão iguais
str(gps)       # Verifique a estrutura de gps
str(PCs)       # Verifique a estrutura de PCs
all(gps %in% rownames(PCs))  # Comparar gps com os nomes das linhas (espécies) de PCs
gp <- factor(gps)  # Assegure-se de que 'gps' está corretamente convertido em fator

gps <- factor(gps)  # Tornar 'gps' um fator
# Reorganizando 'gps' de acordo com os rownames de 'PCs'
gps <- gps[match(rownames(PCs), gps)]
# Confirmando que as espécies em gps estão alinhadas com rownames(PCs)
all(gps == rownames(PCs))  # Isso deve retornar TRUE se os nomes estiverem corretamente alinhados
rownames(PCs) <- gps
# Verificando se os nomes de gps estão nos rownames de PCs
setdiff(gps, rownames(PCs))  # Isso deve retornar NULL ou um vetor vazio


result <- compare.multi.evol.rates(PCs[rownames(PCs)], gps, bats.tree$tip.label, iter = 999)


PCs_filtrado <- PCs[rownames(PCs) %in% bats.tree$tip.label, , drop = FALSE]
gps_filtrado <- gps[rownames(PCs_filtrado)]
all(rownames(PCs_filtrado) %in% bats.tree$tip.label)  # Deve retornar TRUE
result <- compare.multi.evol.rates(PCs[row.names(PCs)], gp = gps_filtrado, phy = bats.tree, iter = 999)
