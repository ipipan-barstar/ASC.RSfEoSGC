
skiplowcardclusters = 0;
extraclusters = 0;
#extraclusters=2; 
#extraclusters=100; 

clToRemove = c();

# is the grqaphics to be saved 
IsPNG = FALSE;
# is the heading of the graphics to be shown
IsMAIN = TRUE
SLAB = 1.5    # font size for labels
SAXIS = 1.6    # font size for axis text
SMAIN = 1        # font size for graphicvs header
SSUB = 0.9    # font size for subtitkles


MAIN <- function(text)
{
    if (IsMAIN)
      return(text)
    else
      return("")
}


while (sink.number() > 0)
  sink();
while (!is.null(dev.list()))
  dev.off();
source("wersje_spectral.R");

cat("

=================================
R_NNN_spectral_clusterings
=================================
Prior usage of R_xxx_read and R_xxx_doc_matrix assumed

input:
nnn.type
  nnn_cats ; 
 RES0; 

output 
ClsTXT
RES

");

if (length(clToRemove) > 0)
  cat("We remove classes ", paste(clToRemove, collapse = " and "), ".");
flush.console();

if (IsPNG)
{
    if (IsMAIN)
      it = ""
    else
      it = "NT"
  combPath = paste0(TARGET_DIR_PREFIX, "Wyniki_", substr(nnn.type, 1, 3), "/", it, nnn.type, "classesComb_")
  normPath = paste0(TARGET_DIR_PREFIX, "Wyniki_", substr(nnn.type, 1, 3), "/", it, nnn.type, "classesNorm_")
  basepath = paste0(TARGET_DIR_PREFIX, "Wyniki_", substr(nnn.type, 1, 3))
  ifelse(!dir.exists(file.path(getwd(), basepath)), dir.create(getwd(), basepath), FALSE)
}


ClsTXT = nnn_cats;
RES = RES0;
RES = RES[!(ClsTXT %in% clToRemove),]
ClsTXT = ClsTXT[!(ClsTXT %in% clToRemove)]

{
  cat("
removing empty rows
");
  er = c();
  st = c();
  for (jr in 1:dim(RES)[1])
  {
    if (sum(RES[jr,]) == 0)
      er = c(er, jr)
    else
      st = c(st, jr);
  }
  print(er);
  RES = RES[st,];
  ClsTXT = ClsTXT[st];
  rownames(RES) = c()
}

classnames = (names(table(ClsTXT)))
trueNoCl = length(classnames);
Cls = c();

     for (j in ClsTXT)
       Cls = c(Cls, which(j == classnames));

cat("
Computing similarity matrix
");
flush.console();

noDocs = dim(RES)[1];
wds = colnames(RES)
wds = wds[substr(wds, 1, 1) == 'W'];
noWords = length(wds);

S = matrix(nrow = noDocs, ncol = noDocs);
S[,] = 0;
S = as.matrix(RES[, wds]) %*% t(as.matrix(RES[, wds]))
D = diag(S)
Ds = sqrt(D);
Dv = Ds %*% t(Ds)
S = S / Dv;


cat("
Comparison of results of combinatorial spectral clustering with original data: - natural length rows 
");
cmpClusterings(ClsTXT, CombinatorialSpectral(S, trueNoCl, unitrow = FALSE)$newcls);
readline("0. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of combinatorial spectral clustering with original data: - unit length rows 
");
cmpClusterings(ClsTXT, CombinatorialSpectral(S, trueNoCl, unitrow = TRUE)$newcls);
readline("1. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of combinatorial spectral clustering with original data: - unit length rows and one additional dimension used
");
cmpClusterings(ClsTXT, CombinatorialSpectral(S, trueNoCl + 1, trueNoCl, unitrow = TRUE)$newcls);
readline("2. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of Kamvar spectral clustering with original data:
");
cmpClusterings(ClsTXT, KamvarSpectral(S, trueNoCl)$newcls);
readline("3. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of Kamvar spectral clustering with original data: with additional dimension used 
");
cmpClusterings(ClsTXT, KamvarSpectral(S, trueNoCl + 1, trueNoCl)$newcls);
readline("4. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();


cat("
Comparison of results of normalized spectral clustering with original data:
");
cmpClusterings(ClsTXT, NormalizedSpectral(S, trueNoCl)$newcls);
readline("5. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of Normalized spectral clustering with original data: - unit length rows 
");
cmpClusterings(ClsTXT, NormalizedSpectral(S, trueNoCl, unitrow = TRUE)$newcls);
readline("6. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();

cat("
Comparison of results of Normalized spectral clustering with original data: - unit length rows and one additional dimension used
");
cmpClusterings(ClsTXT, NormalizedSpectral(S, trueNoCl + 1, trueNoCl, unitrow = TRUE)$newcls);

readline("7. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();


cat("
Comparison of results of Normalized spectral clustering with original data: - unit length rows and one additional dimension used PLUS prior SVD analysis
");
r = 295; cmpClusterings(ClsTXT, NormalizedSpectral(similaritySVDapproximation(S, r), trueNoCl + 1, trueNoCl, unitrow = TRUE)$newcls);

readline("8. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();


cat("
Comparison of results of Normalized spectral clustering with original data: - unit length rows  PLUS prior SVD analysis
");
r = 295; cmpClusterings(ClsTXT, NormalizedSpectral(similaritySVDapproximation(S, r), trueNoCl, unitrow = TRUE)$newcls);

readline("9. see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();


nd1o = NormalizedSpectral(S, trueNoCl + 1, trueNoCl, unitrow = TRUE)
sub = sample(dim(S)[1], round(dim(S)[1] * 0.2))
nd1sub = NormalizedSpectral(S[sub, sub], trueNoCl + 1, trueNoCl, unitrow = TRUE)


readline("10.see the results  and PRESS ENTER"); cat("processing...\n"); flush.console();


cat('
#############################
         THE END', nnn.type, '
#############################

');


############################################
############################################
##########                        ##########
########## AD MAIOREM DEI GLORIAM ##########
##########                        ##########
############################################
############################################

