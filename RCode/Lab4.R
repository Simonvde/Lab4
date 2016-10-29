setwd("~/Desktop/Lab4")

library(igraph)

library("VGAM")
library("stats4")

data <- read.table("./data/Catalan_dependency_tree_metrics.txt",
                              header = FALSE)

plot(data[c(0,2)])

plot(data[c(1,3)])

data_n = length(data[1]$V1)


for(i in 1:data_n){
  row <- data[i,]
  n <- row$V1
  k2 <- row$V2
  delta <- row$V3
  if((4-6/n)>k2 || k2>n-1){
    #False observations are due to rounding errors.
    print(c("FALSE",i,"k2",k2,(4-6/n),n-1))
    print(row)
  } 
  if(k2*n/(8*(n-1))+1/2>delta || delta>n-1){
    print(c("FALSE",i,"delta"))
  }
}



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



# output <- matrix(nrow=length(source$language),ncol=5)
# 
# write_AICs <- function(language,file) {
#   degree_sequence = read.table(file, header = FALSE)
#   #cat(language,AICs(degree_sequence$V1),"\n")
#   return(AICs(degree_sequence$V1))
# }
# 
# for (x in 1:nrow(source)) {
#   output[x,] <- write_AICs(source$language[x], source$file[x])
# }
# dimnames(output) <- list(source$language, c("zeta","zeta_2","RT_zeta","geom","poisson"))
# output
# library(xtable)
# xtable(output)



Catalan = read.table("./data/Catalan_dependency_tree_metrics.txt", header = FALSE)
colnames(Catalan) = c("vertices","degree_2nd_moment", "mean_length")
Catalan = Catalan[order(Catalan$vertices), ]
plot(Catalan$vertices, Catalan$mean_length,
     xlab = "vertices", ylab = "mean dependency length")
plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")

mean_Catalan = aggregate(Catalan, list(Catalan$vertices), mean)

plot(mean_Catalan$vertices, mean_Catalan$mean_length,
     xlab = "vertices", ylab = "mean mean dependency length")
plot(log(mean_Catalan$vertices), log(mean_Catalan$mean_length),
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

#For function with no parameters: the random tree model!
tree_metrics <- Catalan
RSS <- sum((tree_metrics$mean_length-(tree_metrics$vertices+1)/3)^2)

n <- length(tree_metrics$vertices)
p <- 0
s <- sqrt(RSS/(n - p))

AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)



plot(log(Catalan$vertices), log(Catalan$mean_length),
     xlab = "log(vertices)", ylab = "log(mean dependency length)")
lines(log(Catalan$vertices), log(fitted(nonlinear_model)), col = "green")