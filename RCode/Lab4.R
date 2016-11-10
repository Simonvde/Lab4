setwd("~/Desktop/CSN-Labs/Lab4")

library(igraph)

library("VGAM")
library("stats4")
library(xtable)
library("data.table")
#install.packages("nls2")
library("nls2")

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
colnames(Czech) = c("vertices","degree_2nd_moment", "mean_length")
Czech = Czech[order(Czech$vertices), ]
Italian <- read.table("./data/Italian_dependency_tree_metrics.txt",
                              header = FALSE)
colnames(Italian) = c("vertices","degree_2nd_moment", "mean_length")

Arabic <- read.table("./data/Arabic_dependency_tree_metrics.txt",
                              header = FALSE)
colnames(Arabic) = c("vertices","degree_2nd_moment", "mean_length")
Arabic = Arabic[order(Arabic$vertices), ]


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


# Heteroscedasticity (useless) ------------------------------------------------------

#plot(Catalan$vertices, Catalan$mean_length)
#plot(aggregate(list(Catalan$vertices), Catalan$mean_length, mean))
#plot(aggregate(resid(nonlinear_model),list(Catalan$vertices),mean))

#plot(aggregate(resid(nonlinear_model),list(Catalan$vertices),mean))
#plot(fitted(nonlinear_model),abs(residuals(nonlinear_model)))

#plot(Catalan$mean_length,resid(nonlinear_model))
#plot(aggregate(.1/resid(nonlinear_model),list(Catalan$mean_length),mean))

#plot(log(Catalan$vertices), log(Catalan$mean_length),
#     xlab = "log(vertices)", ylab = "log(mean dependency length)")
#lines(log(Catalan$vertices), log(fitted(nonlinear_model)), col = "green")


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
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  print(lmFit)
  coefs <- coef(lmFit)
  b <- coefs[1]
  vals <- c()
  if(length(data$vertices)<200){
    grid <- expand.grid(list(b = seq(b*1,b/1, by = -(b/1)),
                             d = seq(-3,5,by=.5)))
    model.grid <- nls2(mean_length~model1p(vertices,b,d),data=data,start=grid,algorithm="brute-force")
    vals <- model.grid$m$getPars()
    vals[1] <- vals[1]+.43
  }
  else{
    vals <-  c(b,-1)
    names(vals) <- mCall[c("b","d")]
  }
  vals
}

SSmodel1p <- selfStart(model1p, model1pInit,c("b","d"))



model2p <- function(predictor,a,b,d){
  a*predictor^b+d
}

model2pInit <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  vals <- c()
  if(length(data$vertices)<200){
    grid <- expand.grid(list(a=seq(a*2,a/2, by = -(a/4)),
                             b = seq(b*2,b/2, by = -(b/4)),
                             d = seq(-5,3,by=.5)))
    model.grid <- nls2(mean_length~model2p(vertices,a,b,d),data=data,start=grid,algorithm="brute-force")
    vals <- model.grid$m$getPars()
  } else{
    vals <-  c(a,b,d)
    names(vals) <- mCall[c("a","b","d")]
  }
  vals
}

SSmodel2p <- selfStart(model2p, model2pInit,c("a", "b","d"))




model3p <- function(predictor,a,c,d){
  a*exp(c*predictor)+d
}
model3pInit <- function(mCall,LHS,data){
  if(length(data$vertices)<200){
    grid <- expand.grid(list(a=seq(-6,3,by=1),
                             c=seq(-.1,.1,by=.01),
                             d = seq(-2,4,by=1)))
    model.grid <- nls2(mean_length~model3p(vertices,a,c,d),data=data,start=grid,algorithm="brute-force")
    return(model.grid$m$getPars())
  } 
  else{
    value <-  c(-2,-0.1,2.5)
    if(length(data[[1]])==25037) value <- c(a,c,-1000)
    names(value) <- mCall[c("a","c","d")]
    return(value)
  }
  
}
SSmodel3p <- selfStart(model3p, model3pInit,c("a","c","d"))




model4 <- function(predictor,a){
  a*log(predictor)
}
model4Init <- function(mCall,LHS,data){
  value <-  c(1)
  names(value) <- mCall[c("a")]
  value
}
SSmodel4 <- selfStart(model4, model4Init,c("a"))



model4p <- function(predictor,a,d){
  a*log(predictor)+d
}
model4pInit <- function(mCall,LHS,data){
  value <-  c(1,1)
  names(value) <- mCall[c("a","d")]
  value
}
SSmodel4p <- selfStart(model4p, model4pInit,c("a","d"))



model5 <- function(predictor,a,b,c){
  a*(predictor^b)*exp(c*predictor)
}
model5Init <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ (log(xy[, "x"])+xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  c <- coefs[3]
  print(c(a,b,c))

  value <-  c(a,b,c)
  names(value) <- mCall[c("a","b","c")]
  value
}
SSmodel5 <- selfStart(model5, model5Init,c("a","b","c"))



# language = read.table(source$file[5], header = FALSE)
#     colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
#     language = language[order(language$vertices), ]
#     language <- aggregate(language,list(language$vertices),mean)
#     
# getInitial(models.abstract[[9]],data=language)
# nls(models.abstract[[9]],data=language,trace=T,control=nls.control(maxiter=5000))
# nls(models.abstract[[9]],data=language,trace=T,start=list(a=1e-6,b=8,c=0.01))
# plot(language$vertices,language$mean_length)
# 
# for (x in c(1:10)) {
#     language = read.table(source$file[x], header = FALSE)
#     colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
#     language = language[order(language$vertices), ]
#     
#     #language <- aggregate(language,list(language$vertices),mean)
#     
#     model <- nls(mean_length~SSmodel5p(vertices,a,b,c,d),data=language,control=list(maxiter=500))
#     plot(language$vertices,language$mean_length,main=source$language[x])
#     lines(language$vertices, fitted(model), col = 'red')
#     Sys.sleep(.5)
#     #plotAllModelsLog(language,source$language[x])
# }


model5p <- function(predictor,a,b,c,d){
  a*(predictor^b)*exp(c*predictor)+d
}
model5pInit <- function(mCall,LHS,data){
  xy <- sortedXyData(mCall[["predictor"]],LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ (log(xy[, "x"])+xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  c <- coefs[3]
  grid <- expand.grid(list(a=seq(a*2,a/2, by = -(a/2)),
                           b = seq(b*2,b/2, by = -(b/2)),
                           c=seq(c*4,c/2, by = -(c/2)),
                           d = seq(-2,2,by=1)))
  model.grid <- nls2(mean_length~model5p(vertices,a,b,c,d),data=data,start=grid,algorithm="brute-force")
  model.grid$m$getPars()
  
  #value <-  c(a,b,c,0)
  #names(value) <- mCall[c("a","b","c","d")]
  #value
}
SSmodel5p <- selfStart(model5p, model5pInit,c("a","b","c","d"))



models.abstract <- c(mean_length~SSmodel1(vertices,b),
            mean_length~SSmodel2(vertices,a,b),
            mean_length~SSmodel3(vertices,a,c),
            mean_length~SSmodel1p(vertices,b,d),
            mean_length~SSmodel2p(vertices,a,b,d),
            mean_length~SSmodel3p(vertices,a,c,d),
            mean_length~SSmodel4(vertices,a),
            mean_length~SSmodel4p(vertices,a,d),
            mean_length~SSmodel5(vertices,a,b,c),
            mean_length~SSmodel5p(vertices,a,b,c,d))
#models <- lapply(models.abstract,function(x) nls(x,data=Catalan))

modelsLanguage <- function(language,name){
  colnames(language) = c("vertices","degree_2nd_moment", "mean_length")
  language = language[order(language$vertices), ]
  #language <- aggregate(language, list(language$vertices), mean)
  models <- lapply(models.abstract,function(x) nls(x,data=language,control=nls.control(maxiter=5000,warnOnly=T,minFactor=1e-7,tol=1e-1)))
  models
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


modelsValuesTypeList <- list()
valuesType <- function(language,type,name){
  models <- modelsLanguage(language,name)
  vals <- sapply(models,values)
  result <- numeric(length(models)+1)
  result[1] <- tree_metrics(language)[type]
  for(i in 1:length(models)){
    result[i+1] <- (vals[[i]][type])
  }
  result
}


write_AICs <- function(file,name,type) {
    language = read.table(file, header = FALSE)
    colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
    language = language[order(language$vertices), ]
    if(type=="diffAIC"){
      result <- valuesType(language,"AIC",name)
      return(result-min(result))
    }
    return(valuesType(language,type))
}

table2 <- function(type){
  output <- matrix(nrow=length(source$language),ncol=length(models.abstract)+1)
  for (x in 1:nrow(source)) {
    output[x,] <- write_AICs(source$file[x],source$language[x],type)
  }
  dimnames(output) <- list(source$language, c("0","1","2","3","1+","2+","3+","4","4+","5","5+"))
  output
}
table2("s")
table2("AIC")
table2("diffAIC")
xtable(output)

table3 <- function(){
  output <- matrix(nrow=length(source$language),ncol=23)
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

plotAllModels <- function(language,name){
  colors <- c("red","green","blue","darkred","darkgreen")
  plot(language$vertices,language$mean_length,main = name)
  cc=0
  models <- modelsLanguage(language)
  for(x in models){
    cc=cc+1
    lines(language$vertices, fitted(x), col = colors[cc])
  }
}

modelNames <- c("Model 1","Model 2","Model 3","Model 1+","Model 2+","Model 3+",
                "Model 4","Model 4+","Model 5","Model 5+")
#Plot best model in Log scale:
modelsLanguageList <- list()
vectLanguageList <- list()
plotAllModelsLog <- function(language,name){
  models <- modelsLanguage(language)
  vect <- valuesType(language,"AIC")
  indexMinimumModel <- match(min(vect),vect)
  minModel <- models[[indexMinimumModel-1]]

  plot(log(language$vertices),log(language$mean_length),
       main = paste(name,modelNames[indexMinimumModel-1]),
       xlab="log(vertices)",ylab="log(mean_length)")
  lines(log(language$vertices), log(fitted(minModel)), col = "red")
}

plotHomoscedasticity <- function(language, name) {
  dt = data.table(language)
  dt = dt[, list(length = .N, sd=sd(mean_length)), by = vertices]
  plot(dt$vertices, dt$sd)
  abline(lm(sd ~ vertices, data = dt),col="red",lwd=1.5)
}

for (x in 1:nrow(source)) {
    language = read.table(source$file[x], header = FALSE)
    colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
    language = language[order(language$vertices), ]
    
    png(file=paste("images//homoscedasticity_",source$language[x],".png", sep=""))
    plotHomoscedasticity(language,source$language[x])
    dev.off()
    
    language <- aggregate(language,list(language$vertices),mean)
    
    png(file=paste("images//bestModel_",source$language[x],".png",sep=""))
    plotAllModelsLog(language,source$language[x])
    dev.off()
}
for (x in 1:nrow(source)) {
    language = read.table(source$file[x], header = FALSE)
    colnames(language) <- c("vertices","degree_2nd_moment", "mean_length")
    language = language[order(language$vertices), ]
    language <- aggregate(language,list(language$vertices),mean)
    plotAllModelsLog(language,source$language[x])
}
