

while (sink.number() > 0)
  sink();
while (!is.null(dev.list()))
  dev.off();

source("kmeansWITHweights.R");
source("cut_criteria.R");

cat("
################################################################
# R_find_NRcut -  Derivation
# B-based Clustering
################################################################
");


Mprint <- function(M, Mname)
{
  cat("\n Matrix ", Mname, "\n");
  print(dim(M))
  diff = 0; for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2])
      if (abs(M[j, k] - M[k, j]) > diff)
        diff = abs(M[j, k] - M[k, j])
  cat("\n assymetry=", diff, "\n");

  flush.console()
}

mkSYMprint <- function(M, Mname)
{
  cat("\nSymmetrize  Matrix ", Mname, "\n");
  print(dim(M))
  for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2])
    {
      s = (M[j, k] + M[k, j]) / 2; M[j, k] = s;
      M[k, j] = s;
    }
  diff = 0;
  for (j in 1:dim(M)[1])
    for (k in 1:dim(M)[2])
      if (abs(M[j, k] - M[k, j]) > diff)
        diff = abs(M[j, k] - M[k, j])
  cat("\n assymetry=", diff, "\n");
  flush.console()
  return(M)
}

source("spectral_versions.R");
cat("\n Saving the data: RES0 and nnn_cats");
targetDir = paste0(TARGET_DIR_PREFIX, "Results_", substr(nnn.type, 1, 3));
target.data = cbind(theGROUP = nnn_cats, RES0)

if (!dir.exists(targetDir))
{
  dir.create(targetDir);
  cat("\nDirectory ", targetDir, "created\n\n");
}

write.csv(target.data, paste0(targetDir, "/", nnn.type, "_DATA_and_cats.csv"));
flush.console();

n = dim(S)[1]
# S diagonal must be 0
for (j in 1:dim(S)[1])
  S[j, j] = 0

D = matrix(nrow = n, ncol = n)
D[,] = 0
for (j in 1:n)
  D[j, j] = sum(S[j,])

# correction for unrelated elements 
toRemove = which(diag(D) < 1e-10)
if (length(toRemove) > 0)
{
  toRetain = which(diag(D) >= 1e-10)
  nnn_cats = nnn_cats[toRetain]
  S = S[toRetain, toRetain]
  D = D[toRetain, toRetain]
  RES0 = RES0[toRetain,]
  cat("\n the following documents were removed:");
  print(toRemove);

  n = dim(S)[1]
}

Mprint(D, "D")
w = diag(D)
plot(sort(w), main = "Diagonal  elements sorted");
readline(paste("Diagonal sum", sum(w)));
plot(sort(w + 1), main = "Diagonal'  elements sorted");
readline(paste("Diagonal sum ", sum(w + 1)));
plot(sort(1 / (w + 1)), main = "Diagonal' inverted (F)  elements sorted");
F = sum(1 / (w + 1))
readline(paste("Diagonal inverted sum", sum(F)));

IV = 1:n
IV[] = 1
I = diag(IV)

Isq = matrix(nrow = n, ncol = n)
Isq[,] = 1
E = Isq - I

D = D + I
Dm1 = diag(diag(D)^-1)
nL = Dm1^0.5 %*% (D - S) %*% Dm1^0.5
Mprint(nL, "normL")


DSQ = E %*% Dm1^2 + Dm1^2 %*% E - 2 * Dm1 %*% (S) %*% Dm1

K = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)
Mprint(K, "K")
K = mkSYMprint(K, "B")


egK = eigen(K)

LambdaV = egK$values
if (length(LambdaV[LambdaV < 0]) < 2)
  LambdaV[LambdaV < 0] = 0

if (IsPNG)
  png(paste0(targetDir, "/", "EigenvaluesnBbased", nnn.type, ".png"));
plot(LambdaV, main = "Eigenvalues B-based", ylab = "eigenvalue")
if (IsPNG)
  dev.off()

mEv = min(LambdaV)
cat("\nMinimal eigenvalue of B", mEv, "\n")
sigma = 0; if (mEv < 0) { sigma = -2 * mEv + 1e-10
  cat("\nB is corrected by adding to squared distance matrix ", sigma, "\n")
  for (j in 1:dim(DSQ)[1])
    for (k in 1:dim(DSQ)[2])
      if (j != k)
        DSQ[j, k] = DSQ[j, k] + sigma
  Mprint(DSQ, "DSQ corrected")
  K = -0.5 * (I - Isq / n) %*% DSQ %*% (I - Isq / n)
  Mprint(K, "K corrected")
  K = mkSYMprint(K, "K")

  egK = eigen(K)
  frLambdaV = LambdaV
  LambdaV = egK$values
  cat("\nare there negative eigenvalues now?");
  print(LambdaV[LambdaV < 0])
  mEv = min(LambdaV)
  cat("\nretrial Minimal eigenvalue of K", mEv, "in all", length(LambdaV[LambdaV < 0]), "\n")
  if (length(LambdaV[LambdaV < 0]) < 2)
    LambdaV[LambdaV < 0] = 0 }
Lambda = diag(LambdaV)
Mprint(Lambda, "Lambda")
V = egK$vectors
VT = t(V)
Mprint(V %*% Lambda %*% VT, "reconstructed K")
Mprint(round(V %*% Lambda %*% VT - K, digit = 10), "deviation from K")

i = 2
l = which(S[i,] > 0.001)[1]
x_i = sqrt(Lambda) %*% VT[, i]
x_l = sqrt(Lambda) %*% VT[, l]
dstil = sum((x_i - x_l)^2)
cat("\nS[", i, ",", l, "]=", S[i, l]);
cat("\nsquared distance in the embedding: ", dstil)
cat("\ntrue squared distance              ", DSQ[i, l])
cat("\nintended squared distance          ", 1 / D[i, i]^2 + 1 / D[l, l]^2 - 2 * S[i, l] / (D[i, i] * D[l, l]))
if (sigma > 0)
{
  cat("\nsigma corrected squared distance          ", 1 / D[i, i]^2 + 1 / D[l, l]^2 - 2 * S[i, l] / (D[i, i] * D[l, l]) + sigma)
}
#supertest 
if (1 == 0)
  for (i in 1:n)
    for (l in 1:n)
    {
      x_i = sqrt(Lambda) %*% VT[, i]
      x_l = sqrt(Lambda) %*% VT[, l]
      dstil = sum((x_i - x_l)^2)
      DSQil = DSQ[i, l]
      trueV = (D[i, i] + D[l, l] - 2 * S[i, l]) / (D[i, i] * D[l, l])
      if (abs(dstil - DSQil) > 1e-10)
        stop(paste("dstil DSQil", i, l)); if (i != l)
        if (abs(dstil - trueV) > 1e-10)
          stop(paste("dstil trueV", i, l));
    }

x_i = sqrt(Lambda) %*% VT[, i]
XT = sqrt(Lambda) %*% VT
sum(abs(x_i - XT[, i]))
X = t(XT)

trueNoCl = length(table(nnn_cats))
cat("\nTrue clustering distribution\n");
print(table(nnn_cats))


Nbased = NormalizedSpectral(S, trueNoCl, unitrow = FALSE)
clNbased = Nbased$newcls
cat("\nN-based clustering distribution\n");
print(table(clNbased))


L_LambdaV = Nbased$embed$values
if (IsPNG)
  png(paste0(targetDir, "/", "EigenvaluesNbased", nnn.type, ".png"));
plot(L_LambdaV, main = "Eigenvalues L-based (traditional)", col = "green3", ylab = "eigenvalue")
if (IsPNG)
  dev.off()
plot(c(L_LambdaV, LambdaV[1]), main = "Eigenvalues B (black) and N(green) -based (traditional)", col = "green3")
points(LambdaV)

dic = names(table(nnn_cats))
n_cats = rep(NA, length(nnn_cats))
for (j in 1:length(n_cats))
  n_cats[j] = which(nnn_cats[j] == dic)[1]

k = max(n_cats)
cat("\n n=", n, ", k=", k, " F=", F, "  F-k=", F - k)
cat("\n NRCut criterion", NRCut(S, n_cats))
cat("\n Qnrcut criterion", (F - k + NRCut(S, n_cats)))
xxx = kmeansWITHweightsQ(X[,], n_cats, theWeights = diag(D))
cat("\n Qncut criterion from sum to cluster centers", xxx$tot.withinss)
flush.console(); yyy = kmeansWITHweightsQvD(X[,], n_cats, theWeights = diag(D))
cat("\n Qncut criterion from sum to other cluster  elements", yyy$tot.withinss)


p = dim(Nbased$embed$vectors)[2]
clNembed = Nbased$embed$vectors[, p - 1:trueNoCl]
experiment$figex = "Nemb"
plotEmbed(clNembed, clNbased, "N embedding")
experiment$figex = "NembTrue"
plotEmbed(clNembed, nnn_cats, "N embedding versus hashtags")


fitX = kmeansWITHweights(X[, 1:(trueNoCl + 1)], trueNoCl, nstart = 5, theWeights = diag(D))
print(table(fitX$cluster))
clBbased = fitX$cluster

cnd = TRUE
clBembed = X[cnd, 1:trueNoCl]
experiment$figex = "Bemb"
plotEmbed(clBembed, clBbased[cnd], "B embedding")
experiment$figex = "BembTrue"
plotEmbed(clBembed, nnn_cats[cnd], "B embedding versus hashtags")

plotEmbed(clBembed, clNbased[cnd], "B embedding versus N-embedding")
table(clNbased[cnd])


cat("Clustering result comparison ", " B ", kmeansWITHweightsQ(X[, 1:(trueNoCl + 1)], clBbased, theWeights = diag(D))$tot.withinss, " N ", kmeansWITHweightsQ(X[, 1:(trueNoCl + 1)], clNbased, theWeights = diag(D))$tot.withinss, "true", kmeansWITHweightsQ(X[, 1:(trueNoCl + 1)], n_cats, theWeights = diag(D))$tot.withinss, "\n");


cat("Clustering result comparison ", " B ", kmeansWITHweightsQ(X[, 1:(3 * trueNoCl + 1)], clBbased, theWeights = diag(D))$tot.withinss, " N ", kmeansWITHweightsQ(X[, 1:(3 * trueNoCl + 1)], clNbased, theWeights = diag(D))$tot.withinss, "true", kmeansWITHweightsQ(X[, 1:(3 * trueNoCl + 1)], n_cats, theWeights = diag(D))$tot.withinss, "\n");


cat("Clustering result comparison ", " B ", kmeansWITHweightsQ(X[, 1:(9 * trueNoCl + 1)], clBbased, theWeights = diag(D))$tot.withinss, " N ", kmeansWITHweightsQ(X[, 1:(9 * trueNoCl + 1)], clNbased, theWeights = diag(D))$tot.withinss, "true", kmeansWITHweightsQ(X[, 1:(9 * trueNoCl + 1)], n_cats, theWeights = diag(D))$tot.withinss, "\n");


cat("Clustering result comparison ", " B ", kmeansWITHweightsQ(X[, 1:(19 * trueNoCl + 1)], clBbased, theWeights = diag(D))$tot.withinss, " N ", kmeansWITHweightsQ(X[, 1:(19 * trueNoCl + 1)], clNbased, theWeights = diag(D))$tot.withinss, "true", kmeansWITHweightsQ(X[, 1:(19 * trueNoCl + 1)], n_cats, theWeights = diag(D))$tot.withinss, "\n");

if (n > 500 + 1)
  cat("Clustering result comparison ", " B ", kmeansWITHweightsQ(X[, 1:(500 + 1)], clBbased, theWeights = diag(D))$tot.withinss, " N ", kmeansWITHweightsQ(X[, 1:(500 + 1)], clNbased, theWeights = diag(D))$tot.withinss, "true", kmeansWITHweightsQ(X[, 1:(500 + 1)], n_cats, theWeights = diag(D))$tot.withinss, "\n");


cat("\nB-based clustering distribution\n"); print(table(clBbased))


sink(paste0(targetDir, "/", "normclusterStats", nnn.type, ".txt"))
cat("\n TRUE NUMBER OF CLUSTERS: ", trueNoCl, "\n");
print(table(nnn_cats))
cat("B-based clusters\n");
print(table(fitX$cluster))
cat("N-based clusters\n");
print(table(clNbased))

cat("\n Does N-based GSC imply  our B-based\n");
cmpClusterings(clBbased, clNbased, theComment = " Does N-based GSC imply  our B-based GSC")
cat("\n Does our B-based method imply N-based GSC  \n");
cmpClusterings(clNbased, clBbased)

cat("\n Is true clustering implied with our B-based method \n");
cmpClusterings(nnn_cats, clBbased, theComment = "Is true clustering implied with our B-based method", theLabel = "trueB")

cat("\n Is true clustering  implied  with N-based GSC\n");
cmpClusterings(nnn_cats, clNbased, theComment = "Is true clustering implied with N-based GSC", , theLabel = "trueN")

topwords <- function(RES0, clustering)
{
  for (j in names(table(clustering)))
    if (sum(clustering == j) > 1)
    {
      RES0j = RES0[clustering == j,]
      wc = c()
      for (k in 1:dim(RES0j)[2])
      {
        w = sum(RES0j[, k])
        wc = c(wc, w)
      }
      cat("\n Cluster: ", j, "\n");
      ix = sort(wc, index.return = TRUE, decreasing = TRUE)
      print(colnames(RES0)[ix$ix[1:10]])
    }
}

cat("\n True clusters top words \n");
topwords(RES0, nnn_cats)
cat("\n N-based clusterings top words \n");
topwords(RES0, clNbased)
cat("\n B-based clusterings top words \n");
topwords(RES0, clBbased)

dic = names(table(nnn_cats))
n_cats = rep(NA, length(nnn_cats))
for (j in 1:length(n_cats))
  n_cats[j] = which(nnn_cats[j] == dic)[1]


cat("\n NRCUt B based ", NRCut(S, clBbased))
cat("\n RNCUt N based ", NRCut(S, clNbased))
cat("\n NRCUt true  ", NRCut(S, n_cats))
cat("\n");

sink();
print(table(clNbased))

{
  prRCut("clNbased", S, clNbased)
  prRCut("clBbased", S, clBbased)
  prRCut("true", S, nnn_cats)

  prNCut("clNbased", S, clNbased)
  prNCut("clBbased", S, clBbased)
  prNCut("true", S, nnn_cats)


  prNRCut("clNbased", S, clNbased)
  prNRCut("clBbased", S, clBbased)
  prNRCut("true", S, nnn_cats)

  clrandom = sample.int(trueNoCl, dim(S)[1], replace = TRUE, prob = as.numeric(table(nnn_cats)))
  prNRCut("random", S, clrandom)
  prNCut("random", S, clrandom)
  prRCut("random", S, clrandom)

  if (n > 700 + 1)
  {
    fitXX = kmeansWITHweights(X[, 1:700], trueNoCl, nstart = 5, theWeights = diag(D))
    clBbasedX = fitXX$cluster
    print(table(clBbasedX))
    prNRCut("clBbasedX", S, clBbasedX)
  }
}

flush.console();
plot(sort(DSQ), main = "squared distances sorted");

fitXX = kmeansWITHweights(X[, 1:(3 * trueNoCl)], trueNoCl, nstart = 3, theWeights = rep(1, n))
print(table(fitXX$cluster))

print(table(clBbased, nnn_cats))
print(table(clNbased, nnn_cats))
print(table(clNbased, clBbased))

##########################
# AD MAIOREM DEI GLORIAM #
##########################

while (sink.number() > 0)
  sink();
