##############################################################################################################
##                                                                                                          ##
##                   Análise Multivariada: Uma abordagem aplicada utilizando o software R                   ## 
##                                                                                                          ##
##                                           Análise Discriminante                                          ##
##                                                                                                          ##
################################################  LEGENDA  ###################################################
##                                                                                                          ##
##  x1: Área                                                                                                ##
##  x2: Perímetro                                                                                           ##
##  x3: Capacidade                                                                                          ##
##  x4: Comprimento do kernel                                                                               ## 
##  x5: Largura do kernel                                                                                   ##
##  x6: Coeficiente de Assimetria                                                                           ##
##  x7: Comprimento do sulco do kernel                                                                      ##
##  x8: Tipos de Trigo(Kama, Rosa e Canadense)                                                              ##
##                                                                                                          ##
##############################################################################################################

boxM <-function(data, grouping)
{
  if (!inherits(data, c("data.frame", "matrix")))
    stop("'data' must be a numeric data.frame or matrix!")
  if (length(grouping) != nrow(data))
    stop("incompatible dimensions!")
  dname <- deparse(substitute(data))
  data <- as.matrix(data)
  grouping <- as.factor(as.character(grouping))
  p <- ncol(data)
  nlev <- nlevels(grouping)
  lev <- levels(grouping)
  dfs <- tapply(grouping, grouping, length) - 1
  mats <- aux <- list()
  for(i in 1:nlev) {
    mats[[i]] <- cov(data[grouping == lev[i], ])
    aux[[i]] <- mats[[i]] * dfs[i]
  }
  names(mats) <- lev
  pooled <- Reduce("+", aux) / sum(dfs)
  logdet <- log(unlist(lapply(mats, det)))
  minus2logM <- sum(dfs) * log(det(pooled)) - sum(logdet * dfs)
  sum1 <- sum(1 / (dfs - 1)) 
  sum2 <- sum(1 / ((dfs - 1)^2)) 
  Co <- (((2 * p^2) + (3 * p) - 1) / (6 * (p + 1) *
                                        (nlev - 1))) * (sum1 - (1 / sum(dfs)))
  X2 <- minus2logM * (1 - Co)
  dfchi <- (choose(p, 2) + p) * (nlev - 1)
  pval <- pchisq(X2, dfchi, lower.tail = FALSE)
  out <- structure(
    list(statistic = c("Chi-Sq (approx.)" = X2),
         parameter = c(df = dfchi),
         p.value = pval,
         cov = mats, pooled = pooled, logDet = logdet,
         data.name = dname,
         method = " Box's M-test for Homogeneity of Covariance Matrices"
    ),
    class = c("htest", "boxM")
  )
  return(out)
}


##################################### Função discriminante Linear #########################################

## Primeiramente vamos carregar os pacotes
library(mda)
library(mvShapiroTest)
library(MASS)
library(klaR)

## Leitura do banco de dados 
dados = read.table("AD_seeds.txt", header=TRUE)

## Uma análise descritiva destes dados:
summary(dados)

## Análise descritiva dentro de cada tipo
by(dados[, -8], dados$x8, summary)

## Gráfico de dispersão das variáveis tomadas duas a duas
plot(dados[, -8], col = unclass(dados$x8))

## Boxplots
boxplot(x1 ~ x8, dados, main = "x1 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Área')
boxplot(x2 ~ x8, dados, main = "x2 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Perímetro')
boxplot(x3 ~ x8, dados, main = "x3 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Capacidade')
boxplot(x4 ~ x8, dados, main = "x4 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Comprimento do kernel')
boxplot(x5 ~ x8, dados, main = "x5 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Largura do kernel')
boxplot(x6 ~ x8, dados, main = "x6 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Coeficiente de Assimetria')
boxplot(x7 ~ x8, dados, main = "x7 ~ tipo de grão", xlab = 'Tipo de Grão', ylab='Comprimento do sulco do kernel')

## Teste de normalidade multivariada
X = as.matrix(dados[,1:7])
mvShapiro.Test(X)
## p-valor = 4.996 x 10^-10  -  rejeita-se H0

## Teste de homocedasticidade
## H0: as matrizes de covariâncias são homogêneas
boxM(dados[,1:7], dados$x8)
## p-value = ~ 0: rejeita-se H0  -  deve-se avaliar a função discriminante quadrática

## Estimando as funções discriminantes lineares
X <- as.matrix(dados[, 1:7])
head(X)

dadosLDA <- lda(X, dados$x8)
dadosLDA

## Para estimar as taxas de erro 
b <- predict(dadosLDA)
b

## Matriz de confusão
mc = table(b$class,dados$x8)
mc

mmc = confusion(b$class,dados$x8)

##Taxa de erro aparente 
TEA = attr(mmc,"error")
TEA

##Probabilidade global de acerto
(pAC_LDA = 1-TEA )

## Validação cruzada
dadosLDAcv <- lda(dados[, -8], dados$x8, CV = TRUE)

mcCV = confusion(dadosLDAcv$class, dados$x8)
mcCV

TEA_CV = attr(mcCV,"error")
TEA_CV

## Matriz de novas observações 
nX <- data.frame(x1 = c(14.33443,  18.33429,  11.87386,  13.45), 
                 x2 = c(14.29429,  16.13571,  13.24786,  14.09), 
                 x3 = c(0.8800700, 0.8835171, 0.8494086,  0.87), 
                 x4 = c(5.508057,  6.148029,  5.229514,   4.90),
                 x5 = c(3.244629,  3.677414,  2.853771,   3.33),
                 x6 = c(2.667403,  3.644800,  4.788400,   3.32),
                 x7 = c(5.087214,  6.020600,  5.116400,  5.006))                                                            
nX

## Usando a função predict() para classificar essas novas observações
dados.pred = predict(dadosLDA, nX)
dados.pred
dados.pred$x

## Alguns gráficos e interpretação
## Escores dos grãos nas funções discriminantes
x11()
plot(dadosLDA, col = unclass(dados$x8),abbrev=TRUE)

## Para verificar onde estão localizados os novos indivíduos graficamente
points(dados.pred$x,pch=16,col="blue")

## Verificar qual indivíduo pertence a qual classe
plot(predict(dadosLDA)$x[, 1:2], type = "n")
text(predict(dadosLDA)$x[, 1:2], rownames(dados), col = unclass(dados$x8))

## verificar o posicionamento das médias
plot(dadosLDA, dimen=1, type="both")

## gráfico que mostra a classificação de cada indivíduo para cada combinação de duas variáveis. Mostra 
## também as fronteiras e a TEA em cada situação
windows()
partimat(dados$x8~.,data=dados,method="lda") 



################################### Função discriminante quadrática #######################################

dadosQDA <- qda(X, dados$x8)
dadosQDA

## Para estimar as taxas de erro 
b <- predict(dadosQDA)
b

## Matriz de confusão
mmc = confusion(b$class,dados$x8)
mmc

##Taxa de erro aparente 
TEA = attr(mmc,"error")
TEA

##Probabilidade global de acerto
(pAC_QDA = 1-TEA )

## Validação cruzada
dadosQDAcv <- qda(dados[, -8], dados$x8, CV = TRUE)

mcCV = confusion(dadosQDAcv$class,dados$x8)
mcCV

TEA_CV = attr(mcCV,"error")
TEA_CV

## Matriz de novas observações 
nX <- data.frame(x1 = c(14.33443,  18.33429,  11.87386,  13.45), 
                 x2 = c(14.29429,  16.13571,  13.24786,  14.09), 
                 x3 = c(0.8800700, 0.8835171, 0.8494086,  0.87), 
                 x4 = c(5.508057,  6.148029,  5.229514,   4.90),
                 x5 = c(3.244629,  3.677414,  2.853771,   3.33),
                 x6 = c(2.667403,  3.644800,  4.788400,   3.32),
                 x7 = c(5.087214,  6.020600,  5.116400,  5.006))                                                              
nX

## Usando a função predict() para classificar essas novas observações
dados.pred = predict(dadosQDA, nX)
dados.pred

## Gráfico que mostra a classificação de cada indivíduo para cada combinação de duas variaveis. Mostra também as 
## fronteiras e a TEA em cada situação
x11()
partimat(dados$x8~.,data=dados,method="qda") 




