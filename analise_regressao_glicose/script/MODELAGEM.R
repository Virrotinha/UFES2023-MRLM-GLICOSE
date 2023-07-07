
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
library(mice)   #
library(performance)

# base --------------------------------------------------------------------

df <- read.csv2(
  "bases/glicose.csv",
  sep = ",",
  header=TRUE) 


# manipulando os dados ----------------------------------------------------

glimpse(df) 


df <- df %>%
  mutate(
    DiabetesPedigreeFunction = as.numeric(DiabetesPedigreeFunction),
    BMI = as.numeric(BMI)
  ) %>% 
  select(-Outcome) %>% 
  mutate(Glucose = ifelse(Glucose==0, NA, Glucose),
         BloodPressure = ifelse(BloodPressure==0, NA, BloodPressure),
         SkinThickness = ifelse(SkinThickness==0, NA, SkinThickness),
         Insulin = ifelse(Insulin==0, NA, Insulin),
         BMI = ifelse(BMI==0, NA, BMI),
         DiabetesPedigreeFunction = ifelse(DiabetesPedigreeFunction==0, NA, DiabetesPedigreeFunction),
         Age = ifelse(Age==0, NA, Age)) 


# tratando missing values -------------------------------------------------

quantidade_NA <- df %>% is.na() %>% table %>% as.data.frame() %>% 
  janitor::clean_names() %>% filter(x == "TRUE") %>%  select(freq) %>% rename(quantidadeNA = freq)

tabela_NA <- 
  cbind(
    apply(apply(df, 2, is.na), 2, sum), 
    round(apply(apply(df, 2, is.na), 2, sum) / apply(df, 2, length) * 100, 2)
  ) %>% 
  as.data.frame() %>% 
  mutate(tt = V1/sum(V1)) %>% 
  rename("Qntd. de NA" = V1, "% de NA" = V2) 


df %>% na.omit() %>% dim()
df <- df %>% drop_na(Glucose)
df %>% dim

df <- df %>% 
  mice::mice(method = "cart") %>% 
  complete()

#  49 porcento da insulina são na, o ideal seria retirar essa coluna

# substitui os 0 por NA e retirei as linhas cuja a variavel de interesse(Glicose) é NA, mas não
# omiti as 

# há varioas maneiras para substituir os valores faltantes, utilizando a bilioteca MICE
# conseguimos realizar essa substituição 

# utilizando a função com o mesmo nome da biblioteca

# df <- df %>% mutate(SkinThickness = log(SkinThickness))

# todos os métodos que esse pacote disponibiliza para 
# MOSTRAR TABELA DE NA
# FALAR O PQ QUE O 0 É NA, pessoa com 0 de glicose? 0 de 
# menos a covariavel pregnancies que iremos fazer essa substituição, pois a participante da pesquisa pode ter sim
# ter engravidado 0 vezes.
# glicose, pressão sanguine, insulina


# analise exploratoria ----------------------------------------------------

df %>% ggpairs()

boxplot(df)

cor(df)

# if (cor(df) > 0 & cor(df) < 1 ) {
#   cor(df) > 0.9
# } else {
#   cor(df) < -0.9
# }
# cor(df)
# array(cor(df))


# multicolinearidade ------------------------------------------------------

lm(Glucose ~ ., data=df) %>% vif()	# verifique se há algum valor acima de 10

autocor <- df %>% select(-Glucose) %>% cor() %>% eigen()
max(autocor$values)/min(autocor$values) > 1000
cor(df)

install.packages()

# não há multicolinearidade


# seleção de variáveis ----------------------------------------------------

# modelo inicial
model_inicial <- lm(Glucose ~ .^2, data=df)
model_inicial %>% summary()
model_inicial_summary <- model_inicial %>% summary()
model_inicial_summary$coefficients %>% 
  as.data.frame() %>% 
  janitor::clean_names()
  
opt_model_backw_p <- backward(model,alpha=0.1) # método passo atrás
opt_model_backw_p %>% summary()
opt_model_forw_p<- forward(model,alpha=0.1) # método passo a frente
opt_model_forw_p %>% summary()
opt_model_step_p<- stepWise(model,alpha.enter=0.1,alpha.remove=0.1) # método passo a passo
opt_model_step_p %>% summary()


opt_model_backw_aic<- stepAIC(model,direction="backward") # método passo atrás
opt_model_backw_aic %>% summary()
opt_model_forw_aic<- stepAIC(model,direction="forward") # método passo atrás
opt_model_forw_aic %>% summary()
opt_model_step_aic<- stepAIC(model,direction="both") # método passo atrás
opt_model_step_aic %>% summary()
AIC(opt_model_backw_aic)
AIC(opt_model_forw_aic)
AIC(opt_model_step_aic)
## stepwise - aic ----------------------------------------------------------

# opt_model_step_aic<- stepAIC(model,direction="both") # método passo atrás
# opt_model_step_aic %>% summary()
# AIC(opt_model_step_aic)

## todas as regressões -----------------------------------------------------
library(pracma)
set.seed(0505)
opt_model_all_reg<- glmulti(Glucose ~ .,
                            data = df,
                            crit = aic,		# aic, aicc, bic, bicc
                            level = 2,		# 1 sem interacoes, 2 com
                            method = "g",		# "d", ou "h", ou "g"
                            family = gaussian,
                            fitfunction = glm,	# tipo de modelo (lm, glm, etc)
                            confsetsize = 100	# guarde os melhores 100
)
opt_model_all_reg %>% summary()
best_model<- opt_model_all_reg@objects[[1]]
best_model %>% summary()

# save(tabela_NA,
#      best_model, 
#      opt_model_all_reg, file = "modelos/modelo_allregression.rda")
# load(file = "modelos/modelo_allregression.rda")
# attach("modelos/modelo_allregression.rda")

opt_model_all_reg@formulas[[1]]
media_aic_100_melhores <- mean(opt_model_all_reg@crits)

best_model$coefficients
best_model$effects
# save(best_model, opt_model_all_reg, file = "modelos/modelo_allregression.rda")

# performance::compare_performance(
#   opt_model_backw_p,
#   opt_model_forw_p,
#   opt_model_step_p,
#   opt_model_backw_aic,
#   opt_model_forw_aic,
#   opt_model_step_aic,
#   best_model)

resultado <- performance::check_normality(best_model)
plot(resultado)

# verificando pontos atipicos ---------------------------------------------

# distancia de cox, diagonal da matriz hat e valores resíduos padronizados

fit <- best_model
n <- nrow(df)
k <- length(fit$coefficients)

corte.hii <- 2*k/n
corte.cook <- qf(0.5,k,n-k)
corte.stu <- 2

rst <- rstudent(fit)
hii <- hatvalues(fit)
dcook <- cooks.distance(fit)

obs <- 1:n

df.fit <- data.frame(obs,rst,hii,dcook)




## grafico residuos estudentizados --------------------------------------

df.fit %>% ggplot(aes(x=obs,y=rst)) + 
  geom_point() + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Resíduo Estudentizado") + 
  theme_bw()


## grafico alavancagem ---------------------------------------------------

df.fit %>% ggplot(aes(x=obs,y=hii,ymin=0,ymax=hii)) + 
  geom_point() + 
  geom_linerange() + 
  geom_hline(yintercept = corte.hii, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Alavancagem") + 
  theme_bw()


## grafico distancia de cook---------------------------------------------

df.fit %>% ggplot(aes(x=obs,y=dcook,ymin=0,ymax=dcook)) + 
  geom_point() + 
  geom_linerange() +
  geom_hline(yintercept = corte.cook, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Distância de Cook") + 
  theme_bw()



# GRÁFICO GERAL DE DIAGNÓSTICO

texto <- paste("observ:",obs,"\n",
               "resid_stud:",round(rst,2),"\n",
               "alavancagem:",round(hii,2),"\n",
               "D_Cook:",round(dcook,2),"\n",
               "corte_Cook:",round(corte.cook,2))

ggplotly(
  df.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
    geom_point(aes(size=dcook)) + 
    xlim(0, max(max(hii),corte.hii)) + 
    ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
    geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
    geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
    theme_bw() +
    theme(legend.position="none") + 
    xlab("Alavancagem") + 
    ylab("Resíduo Estudentizado"),
  tooltip = c("text")
)

df.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
  geom_point(aes(size=dcook)) + 
  xlim(0, max(max(hii),corte.hii)) + 
  ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
  theme_bw() +
  theme(legend.position="none") + 
  xlab("Alavancagem") + 
  ylab("Resíduo Estudentizado")


qntd_outlier <- df.fit %>% filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) %>% 
  filter(dcook > corte.cook) %>% nrow()# tamanho do plot

qntd_outlier
corte.cook

# normalidade -------------------------------------------------------------

## histograma r studentizados ----------------------------------------------

rst %>% data.frame() %>% ggplot(aes(x=rst)) + 
  geom_histogram(aes(y=..density..)) + 
  geom_density(alpha=.1, fill="blue") +
  theme_bw()


## gráfico quantil quantil -------------------------------------------------

rst %>% data.frame() %>% ggplot(aes(sample=rst)) + 
  stat_qq() + 
  stat_qq_line() +
  theme_bw()


## testes de normalidade ---------------------------------------------------

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

## Tabela de resultados

testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,t6$method)
estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))
valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,t6$p.value)
resultados <- cbind(estt, valorp)
rownames(resultados) <- testes
colnames(resultados) <- c("Estatística", "p")
round(resultados, digits = 4)


# transformação boxcox ----------------------------------------------------

transf_boxcox<- fit %>% boxcox(data=df)
lambda<- transf_boxcox$x[which.max(transf_boxcox$y)]
lambda

if (lambda == 0) {
  df2 <- df %>%  mutate(Glucose = log(Glucose))
} else {
  df2 <- df %>%  mutate(Glucose = (((Glucose)^lambda - 1 )/ lambda))
}

# XXXX Refazendo as etapas por conta da normalidade ----------------------------

# multicolinearidade ------------------------------------------------------

lm(Glucose ~ ., data=df2) %>% vif()	# verifique se há algum valor acima de 10

# seleção de variáveis ----------------------------------------------------

# modelo inicial
model<- lm(Glucose ~ .^2, data=df2)
model %>% summary()


## stepwise - aic ----------------------------------------------------------

opt_model_step_aic<- stepAIC(model,direction="both") # método passo atrás
opt_model_step_aic %>% summary()
AIC(opt_model_step_aic)


## todas as regressões -----------------------------------------------------
set.seed(0505)
opt_model_all_reg<- glmulti(Glucose ~ .,
                            data = df2,
                            crit = aic,		# aic, aicc, bic, bicc
                            level = 2,		# 1 sem interacoes, 2 com
                            method = "g",		# "d", ou "h", ou "g"
                            family = gaussian,
                            fitfunction = glm,	# tipo de modelo (lm, glm, etc)
                            confsetsize = 100	# guarde os melhores 100
)
opt_model_all_reg %>% summary()
best_model<- opt_model_all_reg@objects[[1]]
best_model %>% summary()


performance::compare_performance(
  # opt_model_backw_p,
  # opt_model_forw_p,
  # opt_model_step_p,
  # opt_model_backw_aic,
  # opt_model_forw_aic,
  opt_model_step_aic,
  best_model)



# verificando pontos atipicos ---------------------------------------------

# distancia de cox, diagonal da matriz hat e valores resíduos padronizados

fit <- best_model
n <- nrow(df2)
k <- length(fit$coefficients)

corte.hii <- 2*k/n
corte.cook <- qf(0.5,k,n-k)
corte.stu <- 2

rst <- rstudent(fit)
hii <- hatvalues(fit)
dcook <- cooks.distance(fit)

obs <- 1:n

df2.fit <- data.frame(obs,rst,hii,dcook)


# GRÁFICO - RESÍDUOS ESTUDENTIZADOS

df2.fit %>% ggplot(aes(x=obs,y=rst)) + 
  geom_point() + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Resíduo Estudentizado") + 
  theme_bw()


# GRÁFICO - ALAVANCAGEM

df2.fit %>% ggplot(aes(x=obs,y=hii,ymin=0,ymax=hii)) + 
  geom_point() + 
  geom_linerange() + 
  geom_hline(yintercept = corte.hii, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Alavancagem") + 
  theme_bw()



# GRÁFICO - DISTÂNCIA DE COOK

df2.fit %>% ggplot(aes(x=obs,y=dcook,ymin=0,ymax=dcook)) + 
  geom_point() + 
  geom_linerange() +
  geom_hline(yintercept = corte.cook, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Distância de Cook") + 
  theme_bw()



# GRÁFICO GERAL DE DIAGNÓSTICO

texto <- paste("observ:",obs,"\n",
               "resid_stud:",round(rst,2),"\n",
               "alavancagem:",round(hii,2),"\n",
               "D_Cook:",round(dcook,2),"\n",
               "corte_Cook:",round(corte.cook,2))

ggplotly(
  df2.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
    geom_point(aes(size=dcook)) + 
    xlim(0, max(max(hii),corte.hii)) + 
    ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
    geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
    geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
    theme_bw() +
    theme(legend.position="none") + 
    xlab("Alavancagem") + 
    ylab("Resíduo Estudentizado"),
  tooltip = c("text")
)

df2.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
  geom_point(aes(size=dcook)) + 
  xlim(0, max(max(hii),corte.hii)) + 
  ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
  theme_bw() +
  theme(legend.position="none") + 
  xlab("Alavancagem") + 
  ylab("Resíduo Estudentizado")


qntd_outlier <- df2.fit %>% filter(rst > 2 | rst < -2) %>% 
  filter(hii > corte.hii) %>% 
  filter(dcook > corte.cook) %>% nrow()# tamanho do plot

qntd_outlier
corte.cook


# normalidade -------------------------------------------------------------

## histograma r studentizados ----------------------------------------------

rst %>% data.frame() %>% ggplot(aes(x=rst)) + 
  geom_histogram(aes(y=..density..)) + 
  geom_density(alpha=.1, fill="blue") +
  theme_bw()


## gráfico quantil quantil -------------------------------------------------

rst %>% data.frame() %>% ggplot(aes(sample=rst)) + 
  stat_qq() + 
  stat_qq_line() +
  theme_bw()


## testes de normalidade ---------------------------------------------------

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

## Tabela de resultados

testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,t6$method)
estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))
valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,t6$p.value)
resultados <- cbind(estt, valorp)
rownames(resultados) <- testes
colnames(resultados) <- c("Estatística", "p")
round(resultados, digits = 4)
alpha = 0.01
round(resultados, digits = 4) %>% as.data.frame() %>%
  mutate(hipotese = ifelse(p < alpha, "rejeita H0", "não rejeita H0"))

# heteroscedasticidade ----------------------------------------------------


# teste de Breusch-Pagan --------------------------------------------------
bptest(fit)			

# plot 1 
het <- performance::check_heteroscedasticity(best_model)
plot(het)

# plot 2 
par(mfrow=c(2,2))
plot(fit)
# par(mfrow=1)
# plot(fit,3)

# XXXX TRATANDO A HETEROCEDASDICIDADE -------------------------------------
attach(df2)
fit_ajustado$formula
fit_gamlss<- gamlss(Y ~ 1 + X4 + X6 + X4:X2 + X5:X4 + X6:X5 + X7:X2 + X7:X3 + X7:X4, 
                    sigma.formula = ~ ., 
                    data=df_ajustado)

fit.gamlss %>% summary()
bptest(fit_ajustado)			
bptest(fit_gamlss)
plot(fit.gamlss)

# r squared ajusted
Rsq(fit.gamlss, type = c("both")) # CoxSnell and CraggUhler
with(summary(fit_ajustado), 1 - deviance/null.deviance) # McFadden’s
# https://search.r-project.org/CRAN/refmans/DescTools/html/PseudoR2.html

# método de comparação
fit_ajustado %>% AIC()
fit.gamlss %>% AIC()


# pressupostos do modelo GAMLSS -------------------------------------------

fit.gamlss %>% vif


# PONTOS INFLUENTES/DE ALAVANCA/OUTLIERS  ---------------------------------

V <- fit.gamlss$sigma.x%*%fit.gamlss$sigma.coefficients %>% 
  exp() %>% 
  as.vector() %>% 
  diag()

V<- V^2
X<- fit.gamlss$mu.x
H<- X%*%solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)

n<- nrow(df)    		# número de observações
k<- ncol(X)	 		# k=p+1 (número de coeficientes)

corte.hii<- 2*k/n		# corte para elementos da diagonal de H
corte.cook<- qf(0.5,k,n-k)	# corte para Distância de Cook
corte.stu<- 2			# corte para resíduos estudentizados

rst<- fit.gamlss$residuals	# resíduos estudentizados
hii<- diag(H)	 		# valores da diagonal da matriz H
dcook<- rst*(hii/(1-hii))*(1/k)	# distância de Cook

obs<- 1:n

df.fit<- data.frame(obs,rst,hii,dcook)



# GRÁFICO - RESÍDUOS ESTUDENTIZADOS

df.fit %>% ggplot(aes(x=obs,y=rst)) + 
  geom_point() + 
  geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Resíduo Estudentizado") + 
  theme_bw()



# GRÁFICO - ALAVANCAGEM

df.fit %>% ggplot(aes(x=obs,y=hii,ymin=0,ymax=hii)) + 
  geom_point() + 
  geom_linerange() + 
  geom_hline(yintercept = corte.hii, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Alavancagem") + 
  theme_bw()



# GRÁFICO - DISTÂNCIA DE COOK

df.fit %>% ggplot(aes(x=obs,y=dcook,ymin=0,ymax=dcook)) + 
  geom_point() + 
  geom_linerange() +
  geom_hline(yintercept = corte.cook, color="red", linetype="dashed") + 
  xlab("Observação") + 
  ylab("Distância de Cook") + 
  theme_bw()



# GRÁFICO GERAL DE DIAGNÓSTICO

texto<- paste("observ:",obs,"\n",
              "resid_stud:",round(rst,2),"\n",
              "alavancagem:",round(hii,2),"\n",
              "D_Cook:",round(dcook,2),"\n",
              "corte_Cook:",round(corte.cook,2))

ggplotly(
  
  df.fit %>% ggplot(aes(x=hii,y=rst,text=texto)) +
    geom_point(aes(size=dcook)) + 
    xlim(0, max(max(hii),corte.hii)) + 
    ylim(-max(abs(rst),corte.stu), max(abs(rst),corte.stu)) + 
    geom_hline(yintercept = c(-corte.stu, corte.stu), color="red", linetype="dashed") + 
    geom_vline(xintercept = corte.hii, color="red", linetype="dashed") + 
    theme_bw() +
    theme(legend.position="none") + 
    xlab("Alavancagem") + 
    ylab("Resíduo Estudentizado"),
  tooltip = c("text")
  
)

corte.cook


###############
# NORMALIDADE #
###############

# HISTOGRAMA

rst %>% data.frame() %>% ggplot(aes(x=rst)) + 
  geom_histogram(aes(y=..density..)) + 
  geom_density(alpha=.1, fill="blue") +
  theme_bw()



# GRÁFICO QUANTIL-QUANTIL

rst %>% data.frame() %>% ggplot(aes(sample=rst)) + 
  stat_qq() + 
  stat_qq_line() +
  theme_bw()



# TESTE DE NORMALIDADE

t1 <- ks.test(rst,"pnorm")	#KS
t2 <- lillie.test(rst)		# Lilliefors
t3 <- cvm.test(rst)		# Cramér-von Mises
t4 <- shapiro.test(rst)		# Shapiro-Wilk
t5 <- sf.test(rst)		# Shapiro-Francia
t6 <- ad.test(rst)		# Anderson-Darling

# Tabela de resultados
testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,t6$method)
estt <- as.numeric(c(t1$statistic,
                     t2$statistic,
                     t3$statistic,
                     t4$statistic,
                     t5$statistic,
                     t6$statistic))
valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,t6$p.value)
resultados <- cbind(estt, valorp)
rownames(resultados) <- testes
colnames(resultados) <- c("Estatística", "p")
round(resultados, digits = 4)



# AUTOCORRELAÇÃO ----------------------------------------------------------

fit_gamlss$y
fit_gamlss$mu.coefficients
fit_gamlss$sigma.
modelo_gamlss_summary <- fit_gamlss %>% summary
modelo_gamlss_summary
