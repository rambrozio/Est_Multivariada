##############################################################################################################
##                                                                                                          ##
##                   Análise Multivariada: Uma abordagem aplicada utilizando o software R                   ## 
##                                                                                                          ##
##                                    Análise de Componentes Principais                                     ##
##                                                                                                          ##
##############################################################################################################

##Obs: Cada Observação se refere a um prato avaliado em um festival gastronômico 

## Carregando o arquivo de dados
(pratos = read.table("ACP_PratosSB.txt", header = TRUE, row.names = 1))
(X = as.matrix(pratos))

## Análise exploratória
summary(X)
windows()
boxplot(X)

library('car')
windows()
plot(pratos)
scatterplotMatrix(X, diagonal = 'none')
cor(cbind(pratos))
## Podemos observar a presença de correlação e heterogeneidade de variâncias: a análise deve ser feita 
## com matriz de correlações

## Teste de espericidade de Bartlett
## H0: R = I (Não existe correlação suficiente para aplicação da técnica multivariada)
## Ha: R <> I (Existe correlação suficiente para aplicação da técnica multivariada)
library(psych)
n = dim(X)[1]
R = cor(X)
cortest.bartlett(R, n)
## significativo (p-valor ~= 0). Rejeita-se H0

## Análise de Componentes Principais utilizando a matriz de correlações (Mais aconselhável, neste caso)
acp_R = princomp(X, cor = T)

## Verificando a quantidade de componentes que se deve reter na analise
## Gráfico screeplot 
screeplot(acp_R, type="l")
abline(h=1, col="red")
## critério de Kaiser 
(k = sum(eigen(R)$values>1))
## Pelo gráfico screeplot e pelo critério de Kaiser devemos reter apenas uma componente.

## Proporção da variação explicada
summary(acp_R)
## Variação explicada pela primeira componente 0.7198257

## Loadings (cargas)
loadings(acp_R)

## Para aparecer as cargas suprimidas
unclass(loadings(acp_R))

## Escores
(escores = acp_R$scores)

## Matriz de correlações  entre variáveis originais e componentes principais
(Cor_vo_cp = cor(X, escores))

## Gráfico biplot
biplot(acp_R, xlim=c(-0.5,1), xlab="Comp.1(71,98%Var.Expl.)", ylab="Comp.2(14,63%Var.Expl.)")

## Conclusão Final: A primeira componente principal explica aproximadamente 72% da variação total e de acordo com a Matriz  
## de correlações entre variáveis originais e componentes principais os  pesos das variáveis Sabor, Atendimento, Temperatura,
## Apresentação e Higiene são negativamente altos para essa componente, ou seja, quanto maior a nota dessas variáveis, menor
## é o escore da primeira componente. Então, a primeira componente principal pode ser entendida como um índice global da 
## qualidade do prato de acordo com os juízes. 
## Assim, escores mais baixo na primeira componente indica que o índice de qualidade é melhor, ou seja, quanto menor o escore 
## dessa componente, melhor é o prato. De acordo com a tabela de escores obtida nessa análise, os pratos P1, P2, P8 e P10 
## possuem a melhor qualidade enquanto que o prato de pior qualidade é o P11.


