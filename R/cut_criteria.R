
NRCut <- function(S, clustering)
{
  nc = 0;
  # zero check for S
  ds = sum(abs(diag(S)))
  if (ds > 0)
    stop("NRCut - diagonal sum of S over 0");
  n = dim(S)[1];
  IV = rep(1, n)
  I = diag(IV)
  DV = as.vector(S %*% IV)
  DV = DV + 1
  D = diag(DV)
  cnl = names(table(clustering))
  k = length(cnl)
  for (j in cnl)
  {
    els = as.character(clustering) == j;
    if (sum(els) > 0)
    {
      V = sum(DV[els]);
      cts = sum(S[els, !els])
      nc = nc + cts / V;
    }
  }
  return(nc)
}

RCut <- function(S, clustering)
{
  nc = 0;
  # zero check for S
  ds = sum(abs(diag(S)))
  if (ds > 0)
    stop("RCut - diagonal sum of S over 0");
  n = dim(S)[1];
  IV = rep(1, n)
  I = diag(IV)
  DV = as.vector(S %*% IV)
  DV[] = 1
  D = diag(DV)
  cnl = names(table(clustering))
  k = length(cnl)
  for (j in cnl)
  {
    els = as.character(clustering) == j;
    if (sum(els) > 0)
    {
      V = sum(DV[els]);
      cts = sum(S[els, !els])
      nc = nc + cts / V;
    }
  }
  return(nc)
}

NCut <- function(S, clustering)
{
  nc = 0;
  # zero check for S
  ds = sum(abs(diag(S)))
  if (ds > 0)
    stop("NCut - diagonal sum of S over 0");
  n = dim(S)[1];
  IV = rep(1, n)
  I = diag(IV)
  DV = as.vector(S %*% IV)
  D = diag(DV)
  cnl = names(table(clustering))
  k = length(cnl)
  for (j in cnl)
  {
    els = as.character(clustering) == j;
    if (sum(els) > 0)
    {
      V = sum(diag(D)[els]);
      cts = sum(S[els, !els])
      if (V > 0)
        nc = nc + cts / V;
    }
  }
  return(nc)
}

prRCut <- function(text, S, clustering) {
  cat("\n RCut criterion for", text, "=", RCut(S, clustering))
}

prNCut <- function(text, S, clustering) {
  cat("\n NCut criterion for", text, "=", NCut(S, clustering))
}

prNRCut <- function(text, S, clustering) {
  cat("\n NRCut criterion for", text, "=", NRCut(S, clustering))
}


plotEmbed <- function(embed, clustering, title = NULL)
{
  par = 12
  isPDF = FALSE

  while (par != 0)
  {
    d2 = par %% 10;
    d1 = (par - d2) / 10;
    if (isPDF)
      pdf(paste0(FIG_DIR, experiment$name, "_", experiment$figex, "_", par, ".pdf"))

    plot(embed[, d1], embed[, d2], , asp = 1, main = title, xlab = paste("dim", d1), ylab = paste("dim", d2));
    colset = c("blue3", "green3", "brown", "yellow", "cyan3", "gray", "navy")
    cnl = names(table(clustering))
    k = length(cnl)
    csid = 0
    for (j in cnl)
    {
      els = as.character(clustering) == j;
      csid = csid + 1
      if (sum(els) > 0 && csid <= length(colset))
      {
        points(embed[els, d1], embed[els, d2], col = colset[csid])
      }
    }
    for (j in cnl)
    {
      els = as.character(clustering) == j;
      csid = csid + 1
      if (sum(els) > 100 && csid <= length(colset))
      {
        ids = which(els)
        sids = sample(ids, round(length(ids) / 10))
        points(embed[sids, d1], embed[sids, d2], col = colset[csid])
      }
    }
    if (isPDF)
      dev.off()
    isPDF = FALSE

    parC = readline(paste("look at plot (0 stop, -1 save)", title, "and \nswitch dimensions up to ", dim(embed)[2], " 0 stops >>> "))
    if (nchar(parC) < 1)
      parC = "0";
    if (parC != "-1")
      par = as.numeric(parC) else isPDF = TRUE
  }
}
