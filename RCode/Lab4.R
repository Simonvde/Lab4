 setwd("~/Desktop/CSN-Labs/Lab4")

library(igraph)

library("VGAM")
library("stats4")
library(xtable)

 
write_summary <- function(language,file) {
  tree = read.table(file, header = FALSE)
  cat(language,length(tree$V1),mean(tree$V1),sd(tree$V1),mean(tree$V3),sd(tree$V3),"\n")
}

source = read.table("./data/languages.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)
for (x in 1:nrow(source)) {
  write_summary(source$language[x], source$file[x])
}

# Introduction (guide) ----------------------------------------------------

Czech <- read.table("./data/Czech_dependency_tree_metrics.txt",
                              header = FALSE)

plot(data[c(0,2)])

plot(data[c(1,3)])

data_n = length(data[1]$V1)


for(i in 1:data_n){
  row <- data[i,]
  n <- row$V1
  k2 <- row$V2
  delta <- row$V3
  error <- 0.00001
  if((4-6/n)>k2+error || k2>n-1+error){
    #False observations are due to rounding errors.
    print(c("FALSE",i,"k2",k2,(4-6/n),n-1))
    print(row)
  } 
  if(k2*n/(8*(n-1))+1/2>delta+error || delta>n-1+error){
    print(c("FALSE",i,"delta"))
  }
}
 
 

Catalan = read.table("./data/Catalan_dependency_tree_metrics.txt", header = FALSE)
colnames(Catalan) = c("vertices","degree_2nd_moment", "mean_length")
Catalan = Catalan[order(Catalan$vertices), ]
plot(Catalan$vertices, Catalan$mean_length,
     xlab = "vertices", ylab = "mean dependency length")
plot((log(Catalan$vertices)), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")

mean_Catalan = aggregate(Catalan, list(Catalan$vertices), mean)

plot((mean_Catalan$vertices), mean_Catalan$mean_length,
     xlab = "vertices", ylab = "mean mean dependency length")
plot((log(mean_Catalan$vertices)), log(mean_Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean mean dependency length)")

plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "vertices", ylab = "mean dependency length")
lines(log(mean_Catalan$vertices),log(mean_Catalan$mean_length), col = "green")
lines(log(mean_Catalan$vertices),log((mean_Catalan$vertices+1)/3), col = "red")


a_initial = 4
b_initial = 4
nonlinear_model = nls(mean_length~a*vertices^b,data=Catalan,
                        start = list(a = a_initial, b = b_initial), trace = TRUE)
linear_model = lm(log(mean_length)~log(vertices), Catalan)
a_initial = exp(coef(linear_model)[1])
b_initial = coef(linear_model)[2]
nonlinear_model = nls(mean_length~a*vertices^b,data=Catalan,
                      start = list(a = a_initial, b = b_initial), trace = TRUE)

deviance(nonlinear_model)
AIC(nonlinear_model)
sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model))
coef(nonlinear_model)
coef(nonlinear_model)["a"]
coef(nonlinear_model)["b"]


# Heteroscedasticity ------------------------------------------------------
 
plot(Catalan$vertices,resid(nonlinear_model))
plot(aggregate(resid(nonlinear_model),list(Catalan$vertices),mean))

plot(aggregate(resid(nonlinear_model),list(Catalan$vertices),mean))
plot(fitted(nonlinear_model),abs(residuals(nonlinear_model)))

plot(Catalan$mean_length,resid(nonlinear_model))
plot(aggregate(.1/resid(nonlinear_model),list(Catalan$mean_length),mean))

plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")
lines(log(Catalan$vertices), log(fitted(nonlinear_model)), col = "green")
 

# Models ------------------------------------------------------------------
 
 
model1 <- function(predictor,b){
  (predictor/2)^b
}
model1Init <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  coefs <- coef(lmFit)
  b <- coefs[1]
  value <-  c(b)
  names(value) <- mCall[c("b")]
  value
}
SSmodel1 <- selfStart(model1, model1Init,c("b"))
 

model2 <- function(predictor,a,b){
  a*predictor^b
}
model2Init <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  value <-  c(a,b)
  names(value) <- mCall[c("a","b")]
  value
}
SSmodel2 <- selfStart(model2, model2Init,c("a", "b"))
 
 
model3 <- function(predictor,a,c){
  a*exp(c*predictor)
}
model3Init <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ xy[, "x"])
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  c <- coefs[2]
  value <-  c(a,c)
  names(value) <- mCall[c("a","c")]
  value
}
SSmodel3 <- selfStart(model3, model3Init,c("a","c"))



model1p <- function(predictor,b,d){
  (predictor/2)^b+d
}
model1pInit <- function(mCall,LHS,data){
  d <- -1
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]+d))
  print(lmFit)
  coefs <- coef(lmFit)
  b <- coefs[1]
  value <-  c(b,d)
  names(value) <- mCall[c("b","d")]
  value
}
SSmodel1p <- selfStart(model1p, model1pInit,c("b","d"))

 
 
model2p <- function(predictor,a,b,d){
  a*predictor^b+d
}
 
model2pInit <- function(mCall,LHS,data){
  d <- .1
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]+d))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  value <-  c(a,b,d)
  names(value) <- mCall[c("a","b","d")]
  value
}

SSmodel2p <- selfStart(model2p, model2pInit,c("a", "b","d"))
 
models.abstract <- c(mean_length~SSmodel1(vertices,b),
            mean_length~SSmodel2(vertices,a,b),
            mean_length~SSmodel3(vertices,a,c),
            mean_length~SSmodel1p(vertices,b,d),
            mean_length~SSmodel2p(vertices,a,b,d))

models <- lapply(models.abstract,function(x) nls(x,data=Catalan))

modelsLanguage <- function(language){
  colnames(language) = c("vertices","degree_2nd_moment", "mean_length")
  language = language[order(language$vertices), ]
  models <- lapply(models.abstract,function(x) nls(x,data=language))
}

# Generate tables ---------------------------------------------------------
 #For function with no parameters: the random tree model!
tree_metrics <- function(language){
  RSS <- sum((language$mean_length-(language$vertices+1)/3)^2)
  
  n <- length(language$vertices)
  p <- 0
  s <- sqrt(RSS/(n - p))
  
  AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
  
  result <- c(RSS,AIC,s)
  names(result) <- c("dev","AIC","s")
  return (result)
}

 
 values  <- function(nonlinear_model){
  results = c(deviance(nonlinear_model),AIC(nonlinear_model),sqrt(deviance(nonlinear_model)/df.residual(nonlinear_model)),coef(nonlinear_model))
  names(results) <- c(c("dev","AIC","s"),names(nonlinear_model$m$getPars()))
  return (results)
}
 

valuesType <- function(language,type){
  models <- modelsLanguage(language)
  vals <- sapply(models,values)
  result <- numeric(length(models)+1)
  result[1] <- tree_metrics(language)[type]
  for(i in 1:length(models)){
    result[i+1] <- (vals[[i]][type])
  }
  result
}
#valuesType(Catalan,"AIC")
 


table2 <- function(type){
  output <- matrix(nrow=length(source$language),ncol=length(models.abstract)+1)
  write_AICs <- function(file) {
    language = read.table(file, header = FALSE)
    colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
    language = language[order(language$vertices), ]
    if(type=="diffAIC"){
      result <- valuesType(language,"AIC")
      return(result-min(result))
    }
    return(valuesType(language,type))
  }
  
  for (x in 1:nrow(source)) {
    output[x,] <- write_AICs(source$file[x])
  }
  dimnames(output) <- list(source$language, c("0","1","2","3","1+","2+"))
  output
}
table2("s")
table2("AIC")
table2("diffAIC")
xtable(output)

table3 <- function(){
  output <- matrix(nrow=length(source$language),ncol=10)
  parameternames=""
    
  for (x in 1:nrow(source)) {
    language = read.table(source$file[x], header = FALSE)
    models <- modelsLanguage(language)
    
    result  <- c()
    for(model in models){
      result = c(result,model$m$getPars())
    }
    parameternames=names(result)
    output[x,] <- result
  }
  
  dimnames(output) <- list(source$language, parameternames)
  output
}
table3()
 


# Generate plots ----------------------------------------------------------

 
colors <- c("red","green","blue","darkred","darkgreen")
plot(Catalan$vertices,Catalan$mean_length)
cc=0
for(x in models){
  cc=cc+1
  lines(Catalan$vertices, fitted(x), col = colors[cc])
}
 
 