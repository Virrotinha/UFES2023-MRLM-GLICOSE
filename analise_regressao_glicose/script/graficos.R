rm(list=ls(all=TRUE))
options(scipen = 999) 

load(file = "modelos/base_tratada.rda")
load(file = "modelos/analise_inicial.rda")

boxplot(df_inicial)

# SUMMARIO
summary(df_inicial)

# MINIMO
df_inicial %>% 
  select(-X8) %>% 
apply(2, min)

# MAXIMO
df_inicial %>% 
  select(-X8) %>% 
  apply(2, max)

# FREQUENCIAS VERIFICAR SE SÃO INTEIROS OU NAO
df_inicial %>% 
  apply(2, table)


ggY <- df_inicial %>%
  ggplot( 
    aes(x = Y,
        fill = ifelse(Y == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (Y) Glicose
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks = seq(0,200,20))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))

ggX1 <- df_inicial %>%
  ggplot( 
    aes(x = X1
        # fill = ifelse(X1 == 0, "yes", "no")
        )) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X1) Pregnancies/Gestações
       ") +
  ylab("Frequência") +
  # scale_fill_manual( 
  #   values = c(
  #     "yes"="red",
  #     "no"="#3F3F3F"), 
  # guide = FALSE)+
  scale_x_continuous(breaks = seq(0,18,2))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))


ggX2 <- df_inicial %>%
  ggplot( 
    aes(x = X2,
        fill = ifelse(X2 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X2) BloodPressure/Pressão Sanguínea
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks = seq(0,130,35))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))

ggX3 <- df_inicial %>%
  ggplot( 
    aes(x = X3,
        fill = ifelse(X3 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X3) SkinThickness/Espessura da Pele
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks = seq(0,100,15))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))

ggX4 <-  df_inicial %>%
  ggplot( 
    aes(x = X4,
        fill = ifelse(X4 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X4) Insulin/Insulina
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks = seq(0,800,150))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))
ggX5 <- df_inicial %>%
  ggplot( 
    aes(x = X5,
        fill = ifelse(X5 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X5) BMI/IMC
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks=seq(0,70,by=10))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))
  
ggX6 <- df_inicial %>%
  ggplot( 
    aes(x = X6,
        fill = ifelse(X6 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
  (X6) DiabetesPedigreeFunction/
Diabetes em função da ancestralidade
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"),
    guide = FALSE)+
  scale_x_continuous(breaks=seq(0, 2.5, by=0.5))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))

ggX7 <- df_inicial %>%
  ggplot( 
    aes(x = X7,
        fill = ifelse(X7 == 0, "yes", "no"))) +
  geom_histogram( 
    alpha = 1,
    position = 'identity') +
  theme_bw() +
  xlab("
       (X7) Age/Idade
       ") +
  ylab("Frequência") +
  scale_fill_manual( 
    values = c(
      "yes"="red",
      "no"="#3F3F3F"), 
    guide = FALSE)+
  scale_x_continuous(breaks = seq(0,80,10))+
  theme(axis.text= element_text(size = 13, family = "serif"),
        axis.title= element_text(size = 15, family = "serif"))





# # par(mfrow=c(2,1))
#  par(bty="n", mfrow=c(2,1))
#  boxplot(df_inicial$X5 ,
#         col="#69b3a2" ,
#         xlab="(X5)",
#         horizontal=TRUE)
#  boxplot(df_inicial$X6 ,
#         col="#69b3a2" ,
#         xlab="(X6) ",
#         horizontal=TRUE)

save(ggY, 
     ggX1, 
     ggX2, 
     ggX3,
     ggX4,
     ggX5,
     ggX6,
     ggX7, file = "modelos/graficos.rda")
###########################

## Generate some data to plot
# X <- mtcars %>% group_by( cyl ) %>% summarize( mpg = mean(mpg) ) %>% ungroup
# 
# ## Add a column indicating whether the category should be highlighted
# X <- X %>% mutate( ToHighlight = ifelse( cyl == 6, "yes", "no" ) )
#
## Data looks like this
## A tibble: 3 x 3
##    cyl      mpg ToHighlight
##  <dbl>    <dbl>       <chr>
##1     4 26.66364          no
##2     6 19.74286         yes
##3     8 15.10000          no
#
## Plot the data
# ggplot( X, aes( x = cyl, y = mpg, fill = ToHighlight ) ) +
#   geom_bar( stat = "identity" ) +
#   scale_fill_manual( values = c( "yes"="tomato", "no"="gray" ), guide = FALSE )
#  
