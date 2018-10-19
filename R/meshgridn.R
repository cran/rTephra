meshgridn = function(L){
  #' N-Dimensional Meshgrid
  #' Make a list in which all combinations of inputs occur.
  #' @param L list containing any number of vectors
  #' @keywords misc
  #' @export
  #' @examples
  #' meshgridn(list(1:2, 5:6, 7:9))
  out = list()
  n = sapply(L, length)
  w = 1:length(L)
  for(i in w){
    out[[i]] = rep(rep(L[[i]], rep(prod(n[w<i]), n[i])), prod(n[w>i])) # prod(NULL) == 1, so this works.
  }
  return(out)
}
