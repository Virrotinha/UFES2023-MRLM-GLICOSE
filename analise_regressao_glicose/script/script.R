# limpando a memória ------------------------------------------------------

rm(list=ls(all=TRUE))
options(scipen = 999) 

# bibliotecas -------------------------------------------------------------

library(MASS)		# stepwise (aic) e transf. BoxCox
library(mixlm)		# stepwise (valor-p)
library(glmulti)	# all regression
library(tidyverse)	# manipulacao de dados
library(plotly)		# gráficos interativos
library(GGally)		# gráfico - matriz de correlação
library(car)    	# vif - multicolinearidade
library(nortest)  	# normalidade
library(lmtest)		# homocedasticidade e auto-correlação
library(gamlss)		# incorporando heterocedasticidade
library(nlme)		# incorporando auto-correlação
library(mice)   
library(performance)

# base --------------------------------------------------------------------

df_inicial <- read.csv2(
  "bases/glicose.csv",
  sep = ",",
  header=TRUE) %>% 
  mutate(
    BMI = as.numeric(BMI),
    DiabetesPedigreeFunction= as.numeric(DiabetesPedigreeFunction)
  ) %>% 
  relocate(Glucose)

total_obs <- nrow(df_inicial)
total_cov <- ncol(df_inicial)


tabela_nomes_variaveis <-
  as.data.frame(
    cbind(
      "caracteristica" = colnames(df_inicial),
      "variavel" = c("Y",
                     "X1",
                     "X2",
                     "X3",
                     "X4",
                     "X5",
                     "X6",
                     "X7",
                     "X8"))
  )

colnames(df_inicial) <- tabela_nomes_variaveis$variavel


# par(mfrow = c(2,4))
# attach(df)
# hist(Glucose)
# hist(Pregnancies)
# hist(BloodPressure)
# hist(SkinThickness)
# hist(Insulin)
# hist(as.numeric(BMI))
# hist(as.numeric(DiabetesPedigreeFunction) )
# hist(Age)

# ggarrange

# DESTACAR OS VALORES ZERO

## inspeção dos dados ------------------------------------------------------

# glimpse(df)

# summary(df)

## manipulando os dados ----------------------------------------------------

df_NA <- df_inicial %>%
  select(-X8) %>%
  mutate(Y = ifelse(Y==0, NA, Y),
         X2 = ifelse(X2==0, NA, X2),
         X3 = ifelse(X3==0, NA, X3),
         X4 = ifelse(X4==0, NA, X4),
         X5 = ifelse(X5==0, NA, X5),
         X6 = ifelse(X6==0, NA, X6),
         X7 = ifelse(X7==0, NA, X7)) 

## tratando missing values -------------------------------------------------

quantidade_NA <- df_NA %>%
  is.na() %>%
  table %>%
  as.data.frame() %>%
  janitor::clean_names() %>%
  filter(x == "TRUE") %>%
  select(freq) %>%
  rename(quantidadeNA = freq)

tabela_NA <-
  cbind(
    apply(apply(df_NA, 2, is.na), 2, sum),
    round(apply(apply(df_NA, 2, is.na), 2, sum) / apply(df_NA, 2, length) * 100, 2)
  ) %>%
  as.data.frame() %>%
  mutate(tt = round(V1/sum(V1)*100, 2)) 

tabela_NA <- rbind(tabela_NA, total = c(sum(tabela_NA$V1), "-","-"))

total_obs_removendo_todos_NAs <- df_NA %>% na.omit() %>% nrow()
df_Y_sem_NA <- df_NA %>% drop_na(Y)
total_obs_removendo_todos_NAs_Y <- df_NA  %>% nrow()


# base tratada ------------------------------------------------------------

df_tratado <- df_Y_sem_NA %>%
  mice::mice(method = "cart") %>%
  complete()

## analise exploratoria ----------------------------------------------------

# df_tratado %>% ggpairs()
# boxplot(df_tratado)
# hist(df_tratado$Y)

## gráfico de correlação ---------------------------------------------------

# library(corrplot)
# corrplot::corrplot(
#   cor(df_tratado),
#   type = 'lower',
#   method = "color",
#   cl.ratio = 0.2,
#   number.cex = 0.7,
#   addCoef.col ='black',
#   tl.cex = 0.7,
#   tl.col = 'black',
#   col = COL2('PuOr', 10),
#   addgrid.col = "black",
# )


## gráfico todas as covariaveis --------------------------------------------

# par(mfrow = c(2,4))
# attach(df_tratado)
# hist(Y)
# hist(X1)
# hist(X2)
# hist(X3)
# hist(X4)
# hist(X5)
# hist(X6)
# hist(X7)

save(df_tratado,
     file = "modelos/base_tratada.rda")
save(df_inicial,
     df_NA,
     df_Y_sem_NA, 
     tabela_nomes_variaveis,
     tabela_NA, 
     quantidade_NA,
     total_obs_removendo_todos_NAs,
     total_obs_removendo_todos_NAs_Y,
     file = "modelos/analise_inicial.rda" )

# MULTICOLINEARIDADE ------------------------------------------------------

# **Autovalores da matriz de correlação**
#   
#   ```{r Autovalores Matriz de Correlação das Covariáveis}
# autocor <- 
#   df_tratado %>% 
#   select(-Y) %>% 
#   cor() %>% 
#   eigen()
# autovalores_cor <- autocor$values %>% 
#   round(3) %>% 
#   as.data.frame() %>% 
#   t()
# 
# rownames(autovalores_cor) <- "Autovalores"
# autovalores_cor_multicolineariedade <- max(autovalores_cor)/min(autovalores_cor) 
# 
# knitr::kable(autovalores_cor,
#              # row.names = "Autovalores",
#              caption = "Autovalores da Matrix de Correlação das Covariáveis"
# )
# ```
# 
# Caso o resultado seja maior que 1000, , o que também não ocorreu.


## vif --------

vif <- lm(Y ~ ., data=df_tratado) %>%
  car::vif() # verifique se há algum valor acima de 10
tabela_vif <- as.data.frame(vif) %>%
  mutate(vif_maior_q_10 = ifelse(vif > 10, "sim","não"))
colnames(tabela_vif) <- c("VIF", "VIF > 10")


## autovalor ----
autocor <-
  df_tratado %>%
  select(-Y) %>%
  cor() %>%
  eigen()
autovalores_cor <- autocor$values
multicolineariedade <- max(autovalores_cor)/min(autovalores_cor)
multicolineariedade
# verificar se é acima > 1000

# MODELO 1 ----------------------------------------------------------------


rm(list=ls(all=TRUE))
load(file = "modelos/base_tratada.rda")

## modelo inicial com interações 2 a 2 -------------------------------------
model_01 <-
  lm(Y ~ .^2,
     data = df_tratado)

model_01_summary <-
  model_01 %>%
  summary()

## tabela dos coeficientes -------------------------------------------------
tabela_model_01 <-
  model_01_summary$coefficients %>%
  as.data.frame() %>%
  janitor::clean_names()

tabela_model_01 <- tabela_model_01 %>%
  mutate(significancia = ifelse(pr_t < 0.01, 1, 0)) %>%
  round(4) %>%
  mutate(significancia = ifelse(significancia == 1, "sim", "nao"))

rownames(tabela_model_01) <-
  c("(Intercepto)",
    rownames(tabela_model_01)[2:length(rownames(tabela_model_01))])
colnames(tabela_model_01) <-
  c("Estimativa",
    "Desvio Padrão",
    "T",
    "P-Valor",
    "Estatísticamente Significante?")

# MODELO 1 - SELEÇÃO DE VARIAVEIS -----------------------------------------

## seleção: todas as regressões -----------------------------------------------------

set.seed(0505)
model_01_todas_reg <- glmulti::glmulti(
  Y ~ .,
  data = df_tratado,
  crit = aic,		# aic, aicc, bic, bicc
  level = 2,		# 1 sem interacoes, 2 com
  method = "g",		# "d", ou "h", ou "g"
  family = gaussian,
  fitfunction = glm,	# tipo de modelo (lm, glm, etc)
  confsetsize = 100	# guarde os melhores 100
)
model_01_todas_reg_summary <- model_01_todas_reg %>% summary()

## melhor modelo -----------------------------------------------------------
model_01_final <- model_01_todas_reg@objects[[1]]
model_01_final_summary <- model_01_final %>% summary()

## tabela dos coeficientes modelo todas regressoes -----------------------
tabela_model_01_final <-
  model_01_final_summary$coefficients %>%
  as.data.frame() %>%
  janitor::clean_names()

tabela_model_01_final <-
  tabela_model_01_final %>%
  mutate(significancia = ifelse(pr_t < 0.01, 1, 0)) %>%
  round(4) %>%
  mutate(significancia = ifelse(significancia == 1, "sim", "nao"))

rownames(tabela_model_01_final) <-
  c("(Intercepto)",
    rownames(tabela_model_01_final)[2:length(rownames(tabela_model_01_final))])
colnames(tabela_model_01_final) <-
  c("Estimativa",
    "Desvio Padrão",
    "T",
    "P-Valor",
    "Estatísticamente Significante?")

tabela_model_01_final

save(model_01,
     model_01_final,
     model_01_todas_reg,
     tabela_model_01,
     tabela_model_01_final,
     file = "modelos/modelo01.rda")

# MODELO 1- PRESSUPOSTOS --------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_tratada.rda")
load(file = "modelos/modelo01.rda")

# valores extremos --------------------------------------------------------

n <- nrow(df_tratado) # n de observações
k <- length(model_01_final$coefficients) # quantidade de covariaveis selecionadas

corte.hii <- 2*k/n
corte.cook <- qf(0.5,k,n-k)
corte.stu <- 2

rst <- rstudent(model_01_final) # residuos estudentizados
hii <- hatvalues(model_01_final) # matrix hat
dcook <- cooks.distance(model_01_final) # distancia de cook

obs <- 1:n # índice das observações

outliers_model_01 <- data.frame(obs, rst, hii, dcook) # base com os valores para verficiar os pontos
outliers_model_01_cortes<- as.data.frame(rbind(corte.hii, corte.cook, corte.stu))

outliers_model_01_qntd_s_cook <- 
  outliers_model_01 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) 

outliers_model_01_qntd_c_cook <-
  outliers_model_01 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) %>%
  filter(dcook > corte.cook) %>% 
  nrow()# tamanho do plot

outliers_model_01_texto <- paste("observ:",obs,"\n",
               "resid_stud:",round(rst,2),"\n",
               "alavancagem:",round(hii,2),"\n",
               "D_Cook:",round(dcook,2),"\n",
               "corte_Cook:",round(corte.cook,2))


save(outliers_model_01,
     outliers_model_01_cortes,
     outliers_model_01_qntd_s_cook,
     outliers_model_01_qntd_c_cook,
     outliers_model_01_texto,
     file = "modelos/modelo01_outilers.rda"
     )


# normalidade -------------------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_tratada.rda")
load(file = "modelos/modelo01.rda")

rst <- rstudent(model_01_final) # residuos estudentizados

## testes de normalidade ---------------------------------------------------

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

## Tabela de resultados

testes <- c(
  t1$method, 
  t2$method,
  t3$method, 
  t4$method,
  t5$method,
  t6$method)

estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))

valorp <- c(t1$p.value, 
            t2$p.value,
            t3$p.value, 
            t4$p.value, 
            t5$p.value,
            t6$p.value)

normalidade_model_01 <- cbind(estt, valorp)
rownames(normalidade_model_01) <- testes
normalidade_model_01 <- round(normalidade_model_01, digits = 4) %>% 
  as.data.frame() %>%
  mutate(resultado = ifelse(valorp < 0.01, "rejeita H0", "não rejeita H0"))
colnames(normalidade_model_01) <- c("Estatística", "P-Valor","Resultado")

rst_model_01 <- rst


## transformação boxcox ----------------------------------------------------

transf_boxcox <- 
  model_01_final %>% 
  MASS::boxcox(data=df_tratado)

lambda_model_01 <- transf_boxcox$x[which.max(transf_boxcox$y)]

if (lambda_model_01 == 0) {
  df_ajustado <- df_tratado %>%  mutate(Y = log(Y))
} else {
  df_ajustado <- df_tratado %>%  mutate(Y = (((Y)^lambda_model_01 - 1 )/ lambda_model_01))
}

save(rst_model_01,
     normalidade_model_01,
     lambda_model_01,
     file = "modelos/modelo01_normalidade.rda")

save(df_ajustado,
     file = "modelos/base_ajustada.rda")

# MODELO 2 ----------------------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")

## seleção: todas as regressões -----------------------------------------------------
set.seed(0505)
model_02_todas_reg <- glmulti::glmulti(
  Y ~ .,
  data = df_ajustado,
  crit = aic,		# aic, aicc, bic, bicc
  level = 2,		# 1 sem interacoes, 2 com
  method = "g",		# "d", ou "h", ou "g"
  family = gaussian,
  fitfunction = glm,	# tipo de modelo (lm, glm, etc)
  confsetsize = 100	# guarde os melhores 100
)

model_02 <- model_02_todas_reg@objects[[1]]
model_02_summary <- model_02 %>% summary()

## tabela dos coeficientes modelo todas regressoes ajustado --------------
tabela_model_02 <- 
  model_02_summary$coefficients %>% 
  as.data.frame() %>% 
  janitor::clean_names()

tabela_model_02 <- 
  tabela_model_02 %>% 
  mutate(significancia = ifelse(pr_t < 0.01, 1, 0)) %>% 
  round(4) %>% 
  mutate(significancia = ifelse(significancia == 1, "sim", "nao"))

rownames(tabela_model_02) <- 
  c("(Intercepto)",
    rownames(tabela_model_02)[2:length(rownames(tabela_model_02))])
colnames(tabela_model_02) <-
  c("Estimativa",
    "Desvio Padrão",
    "T",
    "P-Valor",
    "Estatísticamente Significante?")

tabela_model_02

save(
  model_02,
  model_02_todas_reg,
  tabela_model_02,
  file = "modelos/modelo02.rda"
)

# MODELO 2- PRESSUPOSTOS --------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")
load(file = "modelos/modelo02.rda")

# valores extremos -------------------------------------------------------

n <- nrow(df_ajustado) # n de observações
k <- length(model_02$coefficients) # quantidade de covariaveis selecionadas

corte.hii <- 2*k/n
corte.cook <- qf(0.5,k,n-k)
corte.stu <- 2

rst <- rstudent(model_02) # residuos estudentizados
hii <- hatvalues(model_02) # matrix hat
dcook <- cooks.distance(model_02) # distancia de cook

obs <- 1:n # índice das observações

outliers_model_02 <- data.frame(obs, rst, hii, dcook) # base com os valores para verficiar os pontos
outliers_model_02_cortes<- as.data.frame(rbind(corte.hii, corte.cook, corte.stu))

outliers_model_02_qntd_s_cook <- 
  outliers_model_02 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) 

outliers_model_02_qntd_c_cook <-
  outliers_model_02 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) %>%
  filter(dcook > corte.cook) %>% 
  nrow()# tamanho do plot

outliers_model_02_texto <- paste("observ:",obs,"\n",
                                 "resid_stud:",round(rst,2),"\n",
                                 "alavancagem:",round(hii,2),"\n",
                                 "D_Cook:",round(dcook,2),"\n",
                                 "corte_Cook:",round(corte.cook,2))


save(outliers_model_02,
     outliers_model_02_cortes,
     outliers_model_02_qntd_s_cook,
     outliers_model_02_qntd_c_cook,
     outliers_model_02_texto,
     file = "modelos/modelo02_outilers.rda"
)


# normalidade -------------------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")
load(file = "modelos/modelo02.rda")

rst <- rstudent(model_02) # residuos estudentizados

## testes de normalidade ---------------------------------------------------

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

## Tabela de resultados

testes <- c(
  t1$method, 
  t2$method,
  t3$method, 
  t4$method,
  t5$method,
  t6$method)

estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))

valorp <- c(t1$p.value, 
            t2$p.value,
            t3$p.value, 
            t4$p.value, 
            t5$p.value,
            t6$p.value)

normalidade_model_02 <- cbind(estt, valorp)
rownames(normalidade_model_02) <- testes
normalidade_model_02 <- round(normalidade_model_02, digits = 4) %>% 
  as.data.frame() %>%
  mutate(resultado = ifelse(valorp < 0.01, "rejeita H0", "não rejeita H0"))
colnames(normalidade_model_02) <- c("Estatística", "P-Valor","Resultado")

rst_model_02 <- rst

save(rst_model_02,
     normalidade_model_02,
     file = "modelos/modelo02_normalidade.rda")


# MODELO 3 ----------------------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/modelo02.rda")
load(file = "modelos/base_ajustada.rda")

teste_bp_model_02 <- lmtest::bptest(model_02)

## incorporando a heterocedasticidade

attach(df_ajustado)
model_02$formula
model_03 <- gamlss(Y ~ 1 + X4 + X6 + X4:X2 + X5:X4 + X6:X2 + X6:X5 + X7:X2 + X7:X4, 
                        sigma.formula = ~ . -X2-X1-X3-X4-X6, 
                        data=df_ajustado)

# 0.05
## tabela coeficientes modelo gamlss ---------------------------------------
tabela_model_03 <- 
  model_03 %>% 
  summary()

tabela_model_03 <- 
  tabela_model_03 %>% 
  as.data.frame() %>% 
  janitor::clean_names()

tabela_model_03 <- 
  tabela_model_03 %>% 
  mutate(significancia = ifelse(pr_t < 0.01, 1, 0)) %>% 
  round(4) %>% 
  mutate(significancia = ifelse(significancia == 1, "sim", "nao"))

rownames(tabela_model_03) <- 
  c("(Intercepto)",
    rownames(tabela_model_03)[2:length(rownames(tabela_model_03))])

colnames(tabela_model_03) <-
  c("Estimativa",
    "Desvio Padrão",
    "T",
    "P-Valor",
    "Estatísticamente Significante?")

tabela_model_03

save(teste_bp_model_02,
     model_03,
     tabela_model_03,
     file = "modelos/modelo03.rda"
     )


# MODELO 3-PRESSUPOSTOS ---------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")
load(file = "modelos/modelo03.rda")


## valores extremos --------------------------------------------------------

V <- model_03$sigma.x%*%model_03$sigma.coefficients %>% 
  exp() %>% 
  as.vector() %>% 
  diag()

V<- V^2
X<- model_03$mu.x
H<- X%*%solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)

n<- nrow(df_ajustado)    		# número de observações
k<- ncol(X)	 		# k=p+1 (número de coeficientes)

corte.hii<- 2*k/n		# corte para elementos da diagonal de H
corte.cook<- qf(0.5,k,n-k)	# corte para Distância de Cook
corte.stu<- 2			# corte para resíduos estudentizados

rst<- model_03$residuals	# resíduos estudentizados
hii<- diag(H)	 		# valores da diagonal da matriz H
dcook<- rst*(hii/(1-hii))*(1/k)	# distância de Cook

obs<- 1:n

outliers_model_03 <- data.frame(obs, rst, hii, dcook) # base com os valores para verficiar os pontos
outliers_model_03_cortes<- as.data.frame(rbind(corte.hii, corte.cook, corte.stu))


outliers_model_03_qntd_s_cook <- 
  outliers_model_03 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) 

outliers_model_03_qntd_c_cook <-
  outliers_model_03 %>% 
  filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) %>%
  filter(dcook > corte.cook) %>% 
  nrow()# tamanho do plot

outliers_model_03_texto <- paste("observ:",obs,"\n",
                                 "resid_stud:",round(rst,2),"\n",
                                 "alavancagem:",round(hii,2),"\n",
                                 "D_Cook:",round(dcook,2),"\n",
                                 "corte_Cook:",round(corte.cook,2))


save(outliers_model_03,
     outliers_model_03_cortes,
     outliers_model_03_qntd_s_cook,
     outliers_model_03_qntd_c_cook,
     outliers_model_03_texto,
     file = "modelos/modelo03_outilers.rda"
)

# normalidade -------------------------------------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")
load(file = "modelos/modelo03.rda")

rst <- model_03$residuals # residuos estudentizados

## testes de normalidade ---------------------------------------------------

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

## Tabela de resultados

testes <- c(
  t1$method, 
  t2$method,
  t3$method, 
  t4$method,
  t5$method,
  t6$method)

estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))

valorp <- c(t1$p.value, 
            t2$p.value,
            t3$p.value, 
            t4$p.value, 
            t5$p.value,
            t6$p.value)

normalidade_model_03 <- cbind(estt, valorp)
rownames(normalidade_model_03) <- testes
normalidade_model_03 <- round(normalidade_model_03, digits = 4) %>% 
  as.data.frame() %>%
  mutate(resultado = ifelse(valorp < 0.01, "rejeita H0", "não rejeita H0"))
colnames(normalidade_model_03) <- c("Estatística", "P-Valor","Resultado")

rst_model_03 <- rst

save(rst_model_03,
     normalidade_model_03,
     file = "modelos/modelo03_normalidade.rda")

# autocorrelação ???????????????? -------------------------------

rm(list=ls(all=TRUE))
load(file = "modelos/base_ajustada.rda")
load(file = "modelos/modelo03.rda")

teste_bg_model_03 <- lmtest::bgtest(model_03)

save(teste_bg_model_03, 
     file = "modelos/modelo03_autocorrelacao.rda")
