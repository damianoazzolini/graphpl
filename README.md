# graphpl
SWI prolog package to manage graphs. 

### Basic info
A graph, in this package is represented in this way: `graph(ListOfVertices,ListOfEdges)` where `ListOfVertices` is a list of integer and `ListOfEdges` is a list of predicates `edge/2` or `edge/3` where `edge/2` is used to for unweighted graphs (`edge(NodeA,NodeB)`) and `edge/3` is used for wheighted graphs (`edge(NodeA,NodeB,Cost)`).

### Available Predicates
* `generate_undirected_unweighted_graph/2`
* `generate_undirected_weighted_graph/2`
* `generate_unweighted_graph/2`
* `generate_weighted_graph/2`
* `find_path_unweighted/4`
* `find_path_weighted/5`
* `generate_kn/2`
* `generate_kn_weighted/4`
* `generate_kn_from_vertices/2`
* `cycle_unweighted/3`
* `cycle_weighted/4`
* `is_connected/1`
* `node_degree/3`
* `node_degree_list/2`
* `empty_unweighted_graph/3`
* `empty_weighted_graph/3`
* `is_graph_node/2`
* `is_isolated_node/2`
* `is_graph_edge/2`
* `get_adjacent_nodes/3`
* `graph_reverse_edges/2`
* `spanning_tree/2`
* `mst_prim/3`

### Example

    :- use_module(library(graph)).
    
    test(G):-
	    generate_kn(4,G).
    
    ?- test(T).
    T = graph([1, 2, 3, 4], [edge(1, 2), edge(1, 3), edge(1, 4), edge(2, 3), edge(2, 4), edge(3, 4)]).
    
### Contribution
Feel free to open an issue if you found some problems or pull request if you want to contribute. Feel free also to suggest predicates that could be good to have.
