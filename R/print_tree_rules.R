get_cluster_mapping <- function(tree) {
    # list of factor vectors, each element is a node and each elem in the
    # factor is the majority prediction
    predictions <- partykit::nodeapply(
      partykit::as.simpleparty(tree),
      ids = nodeids(tree, terminal = TRUE),
      function(x) partykit::info_node(x)$prediction)
    
    # Named vector, each name the node id and each element the predicted cluster
    cluster_mapping <- wrapr::named_map_builder(
      names(predictions),
      unlist(predictions))
    
    return(cluster_mapping)
}

get_concensus_rules <-  function(tree) {
  
  cluster_mapping <- get_cluster_mapping(tree)
  
  # List of lists of length 1, each sub-list being a named numeric vector
  number_per_node <- partykit::nodeapply(
    partykit::as.simpleparty(tree),
    ids = nodeids(tree, terminal = TRUE),
    function(x) partykit::info_node(x)$distribution)
  
  # Named vector, each name the node id and each element number of elements for
  # that predicted cluster
  number_per_node <- purrr::map2_dbl(
    number_per_node, 
    cluster_mapping,
    .f = function(x,y) {
      x[[y]][[1]]} )
  
  # Character vector of the rules
  tree_rules <- partykit:::.list.rules.party(tree)
  

  # This is a list whose names are the cluster names and the elements the 
  # name/number of the nodes
  nodes_per_clus <- lapply(split(cluster_mapping, cluster_mapping), names)
  
  # This returns the (overwhelming) rules for each cluster
  rules_per_clus <- lapply(nodes_per_clus, function(x) tree_rules[x])
  
  
  get_majority_rules <- function(rules_in_cluster) {
    split_rules <- strsplit(rules_in_cluster, "\\s&\\s")
    
    all_rules <- unique(unlist(split_rules))
    
    count_per_rule <- function(nodes_logical, number_per_node) {
      # > nodes_logical
      #   19   20   21   24   25   26 
      # TRUE TRUE TRUE TRUE TRUE TRUE 
      # > number_per_node
      #  19  20  21  24  25  26 
      #  59 571  15  48   7   5 
      sum(number_per_node[names(nodes_logical)[which(nodes_logical)]])
      }
    
    .in <- function(x, y) {
      x %in% y
    }
    
    elems_per_rule <- purrr::map(
      all_rules, 
      function(x) purrr::map_lgl(split_rules, function(y) { x %in% y })) %>%
      purrr::map_dbl(count_per_rule, number_per_node)
    
    rules_in_majority <- all_rules[
      elems_per_rule > (max(elems_per_rule)/2) &
        elems_per_rule != max(elems_per_rule)]
    rules_in_all <- all_rules[elems_per_rule == max(elems_per_rule)]
    
    list(all = rules_in_all, majority = rules_in_majority)  
  }
  
  
  purrr::map(rules_per_clus, get_majority_rules)
  # returns a list[each cluster][majority or ]
}
