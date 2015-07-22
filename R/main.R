#' Statistical tests to compare word clouds
#'
#' \code{\link{wordcloudstats}} is the main user interface to conduct a quantitative analysis of
#' the significance of the observed differences between the two word clouds.
#'
#' @import wordcloud
#' 
#' @name wordcloudanalysis-package
#' @docType package
NULL

#' @export
#' @title Quantitative statistical analysis of differences between two word clouds.
#' @description Main user interface function for quantitative statistical analysis between two
#' word clouds. If \code{avoid.words} is specified, then words that ought to be ignored are
#' excluded from the analysis. The common terms between \code{text1} and \code{text2} are
#' considered for statistical analysis. The analysis computes the significance of the
#' observed differences in term occurrence in the \code{text1} and \code{text2}.
#'
#' @param text1 A character vector representing terms corresponding to one word cloud
#' @param text2 A character vector representing terms corresponding to the other word cloud
#' @param names A character vector of length 2 with the two names for the provided lists. The
#'              first element corresponds to \code{text1} and the second element corresponds to \code{text2}.
#' @param avoid.words A character vector of words to ignore when analyzing the two word clouds.
#'                    Defaults to \code{NULL}.
#'
#' @return A list of three elements. \code{outputs} contains a matrix of the significance and
#' counts of each of the common terms between the lists. Specifically, the columns consist of:
#' \itemize{
#' \item count.in.\code{names[1]} (denoting the count of the term in group 1)
#' \item unmatched.term.count.\code{names[1]} (denoting \code{length(text1)} - \code{count.in.names[1]}, or the number
#' of terms in group 1 that are not the specific term)
#' \item count.in.\code{names[2]} (denoting the count of the term in group 2)
#' \item unmatched.term.count.\code{names[2]} (denoting \code{length(text2)} - \code{count.in.names[2]}, or the number
#' of terms in group 2 that are not the specific term)
#' \item p.value (the nominal p-value of the observed differences in the specific term)
#' \item bh.value (the Benjamini-Hochberg adjusted p-value)
#' \item greater.frequency.in (the group in which the frequency of the term is greater: \code{names[1]}, \code{names[2]}, or 'equal')
#' } 
#' 
#' \code{counts} contains a table of
#' the counts of each term. \code{frequency} contains a table of the frequencies of each term
#' in each list.
#'
#' @examples
#' group1 <- c('head','toe','hand',rep('knee',4))
#' group2 <- c(rep('toe',3),'hand',rep('head',2))
#'
#' stats <- wordcloudstats(group1, group2, names=c('Group 1','Group 2'))
#'
#' head(stats$outputs)
#'
#' @seealso \code{\link[stats]{prop.test}}
wordcloudstats <- function(text1, text2, names, avoid.words=NULL) {

  # Remove words that the user wants to avoid
  if (! is.null(avoid.words)) {
    text1 <- text1[!text1 %in% avoid.words]
    text2 <- text2[!text2 %in% avoid.words]
  }

  # Keep only the terms that appear at least once in each list
  common <- intersect(unique(text1),unique(text2))
  text1 <- text1[text1%in%common]
  text2 <- text2[text2%in%common]

  # Basic frequency table
  text1.freq <- unlist(table(text1))/length(text1)
  text2.freq <- unlist(table(text2))/length(text2)

  freq <- matrix(NA,nrow=length(text1.freq),ncol=2)
  rownames(freq) <- names(text1.freq)
  freq[,1] <- text1.freq
  freq[names(text2.freq),2] <- text2.freq
  colnames(freq) <- names

  # Basic count table
  text1.ct <- unlist(table(text1))
  text2.ct <- unlist(table(text2))

  count <- matrix(NA,nrow=length(text1.ct),ncol=2)
  rownames(count) <- names(text1.ct)
  count[,1] <- text1.ct
  count[names(text2.ct),2] <- text2.ct
  colnames(count) <- names

  pvals <- rep(1,nrow(count))
  names(pvals) <- rownames(count)
  g.in <- rep(1,nrow(count))
  names(g.in) <- rownames(count)

  text1.in <- rep(1,nrow(count))
  names(text1.in) <- rownames(count)
  text1.out <- rep(1,nrow(count))
  names(text1.out) <- rownames(count)

  text2.in <- rep(1,nrow(count))
  names(text2.in) <- rownames(count)
  text2.out <- rep(1,nrow(count))
  names(text2.out) <- rownames(count)

  # Run proportions test
  for (term in rownames(count)) {
    text1.yes <- count[term,1]
    text2.yes <- count[term,2]
    text1.no <- length(text1) - text1.yes
    text2.no <- length(text2) - text2.yes

    text1.in[term] <- text1.yes
    text1.out[term] <- text1.no
    text2.in[term] <- text2.yes
    text2.out[term] <- text2.no

    mat <- matrix(NA, nrow=2, ncol=2)
    mat[1,1] <- text1.yes
    mat[1,2] <- text1.no
    mat[2,1] <- text2.yes
    mat[2,2] <- text2.no

    ptest <- prop.test(mat)
    pvals[term] <- ptest$p.value
    maxval<- max(as.numeric(unlist(ptest$estimate)))
    g.in.ind <- grep(maxval,as.numeric(unlist(ptest$estimate)))
    if (length(g.in.ind) > 1)
      g.in[term] <- 'equal'
    else
      g.in[term] <- names[g.in.ind[1]]
  }

  p.adj <- p.adjust(unname(pvals),method='fdr')
  mat.val <- cbind(text1.in,text1.out,text2.in,text2.out,pvals,p.adj,g.in)
  present.names <- paste('count.in',names,sep='.')
  absent.names <- paste('unmatched.term.count',names,sep='.')
  mat.names <- c(present.names[1],absent.names[1], present.names[2],absent.names[2],'p.Value','bh.Value','greater.frequency.in')
  colnames(mat.val) <- mat.names

  list(outputs=mat.val,counts=count,frequency=freq)
}

#' @export
#' @title Qualitative plot of comparison for two word clouds
#' @description Main user interface for qualitative analysis of two word clouds.
#' Creates a plot with frequencies of each group on each axis. The size of each
#' term is inversely related to the q-values provided. The color of each term
#' indicates which group had the greater number of counts for the term occurrence.
#' 
#' Note that all parameters must describe the same number of terms.
#' 
#' @param freq A numeric matrix with two columns. The first column indicates the
#'             frequency of a term in the first group, and the second column indicates
#'             frequency of a term in the second group. The rows describe terms. Rownames
#'             must be set to the terms.
#' @param counts A numeric matrix in the same format as \code{freq} except the entries
#'               represent count values for each of the groups.
#' @param qvals A vector of q-values for each of the terms. The names for the vector
#'              must be the terms. Q-values of 0 are rounded to 0.00001. Note that there
#'              must be at least two distinct q-values (i.e., they cannot all be 1).
#' @param direction A vector of numeric values indicating whether the term belongs to
#'                  group 1 or 2. The names of the vector must be the terms.
#' @param colors A character vector of length three. The first element represents the color
#'               that the term should be plotted in if the count value is greater in the 
#'               first group compared to the second group. The second element represents
#'               the color for terms with greater counts in the second group compared to
#'               the first group. The third element is the color for terms that have equivalent
#'               counts in both groups.Defaults to \code{c('red','blue','black')}
#' @param size.limits Bounds for the sizes of the words. Must be provided in the format max:min.
#'                    The sizes are interpreted by \code{cex}, so they are scaled. Defaults to \code{50:1/20}.
#' @param xlim A numeric vector of length two with the lower and upper bounds for the x-axis of the graph.
#'             Defaults to \code{c(-0.5,1)}.
#' @param ylim A numeric vector of length two with the lower and upper bounds for the y-axis of the graph.
#'             Defaults to \code{c(-0.5,1)}.
#'                    
#' @examples
#' group1 <- c('head','toe','hand',rep('knee',4))
#' group2 <- c(rep('toe',3),'hand',rep('head',2))
#'
#' stats <- wordcloudstats(group1, group2, names=c('Group 1','Group 2'))
#' qval <- runif(nrow(stats$frequency),0.0,1.0)
#' names(qval) <- rownames(stats$outputs)
#' comparisonplot(stats$frequency, stats$counts, qval, colors=c('blue','green','red'))
#'                    
#' @seealso \code{\link[wordcloud]{wordlayout}}
comparisonplot <- function(freq, counts, qvals, colors=c('red','blue','black'), size.limits=50:1/20,
                           xlim=c(-0.5,1),ylim=c(-0.5,1)) {
  
  # Transform q-values
  terms <- names(qvals)
  qvals[qvals==0] <- 0.00001
  qvals <- -log10(qvals)+1
  names(qvals) <- terms
  
  # Ensure everything is in the same order by terms
  qvals <- qvals[rownames(freq)]
  counts <- counts[rownames(freq),]
  
  # Scale q-values to size.limits
  message('Scaling q-values...')
  mn.sc <- min(size.limits)
  mx.sc <- max(size.limits)
  mn.q <- min(qvals)
  mx.q <- max(qvals)
  scaled.vals <- abs(unlist(lapply(qvals,function(x){
    ((x-mn.q) * (mx.sc-mn.sc)) / (mx.q-mn.q) + mn.sc
  })))
  
  # Generate colors
  message('Computing colors...')
  applied.colors <- apply(counts,1,function(x){
    if(x[1] > x[2])
      return(colors[1])
    else if (x[2] > x[1])
      return(colors[2])
    else
      return(colors[3])
  })
  
  # Generate plot
  message('Plotting results...')
  plot(freq[,1],freq[,2],type='n',xlim=xlim,ylim=ylim)
  wl <- wordcloud::wordlayout(freq[,1],freq[,2],rownames(freq),cex=scaled.vals)
  text(wl[,1]+.5*wl[,3], wl[,2]+.5*wl[,4], rownames(freq), cex=scaled.vals, col=applied.colors)
}
