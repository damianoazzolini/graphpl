% graph term form, default
%  Graph = graph([a,b,c,d,e,f,g,h], [e(a,b), e(b,a), … ]).

%:- module(graphpl, [
%    make_undirected_unweighted_graph/2,
%    make_undirected_weighted_graph/2,
%    make_unweighted_graph/2,
%    make_weighted_graph/2,
%    find_path_unweighted/4
%    find_path_weighted/5,
%    make_kn/2,
%    make_kn_weighted/4,
%    make_kn_from_vertices/2,
%    cycle_unweighted/3,
%    cycle_weighted/4,
%    is_connected/1,
%    node_degree/3,
%    node_degree_list/2,
%    empty_unweighted_graph/3,
%    empty_weighted_graph/3,
%    is_graph_node/2,    
%    is_isolated_node/2,
%    is_graph_edge/2,
%    get_adjacent_nodes/3,
%    graph_reverse_edges/2,
%    spanning_tree/2,
%    mst_prim/3
%]).

% make_undirected_unweighted_graph(+ListOfEdges,-Graph) creates a graph in graph-term form
% also removes duplicates edge(a,b) means also edge(b,a)
% undirected and unweighted graph
make_undirected_unweighted_graph(ListOfEdges,graph(VSorted,ListOfEdgesSorted)):-
    findall(X,(member(edge(X,_),ListOfEdges)),VX),
    findall(Y,(member(edge(_,Y),ListOfEdges)),VY),
    append(VX,VY,V),
    sort(V,VSorted),
    remove_duplicated_edges(ListOfEdges,ListOfEdgesNoDup),
    sort(ListOfEdgesNoDup,ListOfEdgesSorted).

remove_duplicated_edges([],[]).
remove_duplicated_edges([edge(X,Y)|T],L):-
    member(edge(Y,X),T),!,
    remove_duplicated_edges(T,L).
remove_duplicated_edges([edge(X,Y)|T],[edge(X,Y)|TE]):-
    remove_duplicated_edges(T,TE).


% make_undirected_weighted_graph(+ListOfEdges,-Graph) creates a graph in graph-term form
% also removes duplicates edge(a,b) means also edge(b,a)
% undirected and weighted graph
make_undirected_weighted_graph(ListOfEdges,graph(VSorted,ListOfEdgesSorted)):-
    findall(X,(member(edge(X,_,_),ListOfEdges)),VX),
    findall(Y,(member(edge(_,Y,_),ListOfEdges)),VY),
    append(VX,VY,V),
    sort(V,VSorted),
    remove_duplicated_edges_w(ListOfEdges,ListOfEdgesNoDup),
    sort(ListOfEdgesNoDup,ListOfEdgesSorted).

remove_duplicated_edges_w([],[]).
remove_duplicated_edges_w([edge(X,Y,V)|T],L):-
    member(edge(Y,X,V),T),!,
    remove_duplicated_edges_w(T,L).
remove_duplicated_edges_w([edge(X,Y,V)|T],[edge(X,Y,V)|TE]):-
    remove_duplicated_edges_w(T,TE).


% make_unweighted_graph(+ListOfEdges,-Graph) creates a graph in graph-term form
% directed and unweighted graph
make_unweighted_graph(ListOfEdges,graph(VSorted,ListOfEdges)):-
    findall(X,(member(edge(X,_),ListOfEdges)),VX),
    findall(Y,(member(edge(_,Y),ListOfEdges)),VY),
    append(VX,VY,V),
    sort(V,VSorted).


% make_weighted_graph(+ListOfEdges,-Graph) creates a graph in graph-term form
% directed and weighted graph
make_weighted_graph(ListOfEdges,graph(VSorted,ListOfEdges)):-
    findall(X,(member(edge(X,_,_),ListOfEdges)),VX),
    findall(Y,(member(edge(_,Y,_),ListOfEdges)),VY),
    append(VX,VY,V),
    sort(V,VSorted).


% find_path_unweighted(+Graph,+V1,+V2,-Path)
% find all path between V1 and V2 in directed unweighted graph
find_path_unweighted(graph(Vertices,Edges),V1,V2,Path):-
    member(V1,Vertices),
    member(V2,Vertices),
    find_path_unweighted_(V1,V2,Edges,[V1],Path).

find_path_unweighted_(V,V,_,L,L).
find_path_unweighted_(V1,V2,Edges,L,Path):-
    % connected(V1,V2,Edges),
    member(edge(V1,V2),Edges),    
    append(L,[V2],Path).
find_path_unweighted_(V1,V2,Edges,L,Path):-
    % connected(V1,VX,Edges),
    member(edge(V1,VX),Edges),
    VX \== V2,
    \+member(VX,L),
    append(L,[VX],L1),
    find_path_unweighted_(VX,V2,Edges,L1,Path).

connected(V1,V2,Edges):-
    member(edge(V1,V2),Edges);
    member(edge(V2,V1),Edges).

connected_w(V1,V2,Edges,C):-
    member(edge(V1,V2,C),Edges);
    member(edge(V2,V1,C),Edges).

% find_path_weighted(+Graph,+V1,+V2,-Path,-Cost)
% find all path between V1 and V2 in a weighted directed graph
find_path_weighted(graph(Vertices,Edges),V1,V2,Path,Cost):-
    member(V1,Vertices),
    member(V2,Vertices),
    find_path_weighted_(V1,V2,Edges,[V1],Path,0,Cost).

find_path_weighted_(V,V,_,L,L,C,C).
find_path_weighted_(V1,V2,Edges,L,Path,CurrC,TotalC):-
    % connected_w(V1,V2,Edges,Cost),
    member(edge(V1,V2,Cost),Edges),    
    append(L,[V2],Path),
    TotalC is CurrC + Cost.
find_path_weighted_(V1,V2,Edges,L,Path,CurrC,TotalC):-
    % connected_w(V1,VX,Edges,Cost),
    member(edge(V1,VX,Cost),Edges),        
    VX \== V2,
    \+member(VX,L),
    append(L,[VX],L1),
    Curr1 is CurrC + Cost,
    find_path_weighted_(VX,V2,Edges,L1,Path,Curr1,TotalC).


% make_kn(+Size,-Graph)
% creates a Kn Graph of size Size with verteces name 1,2,..,N
% only undirected and unweighted graph
make_kn(Size,graph(LV,Comb)):-
    Size1 is Size+1,
    make_ordered_list(1,Size1,LV),
    find_all_combinations(LV,[],Comb).

make_ordered_list(Max,Max,[]):- !.
make_ordered_list(I,Max,[I|T]):-
    I < Max,
    I1 is I+1,
    make_ordered_list(I1,Max,T).

find_all_combinations([_],C,C):- !.
find_all_combinations([H|T],CT,CO):-
    find_combinations(H,T,C),
    append(CT,C,C1),
    find_all_combinations(T,C1,CO).

find_combinations(_,[],[]):- !.
find_combinations(E,[H|T],[edge(E,H)|TE]):-
    find_combinations(E,T,TE).


% make_kn_weighted(+Size,+MinValue,+MaxValue,-Graph)
% creates a Kn undirected Graph of size Size with verteces name 1,2,..,N
% and assign each cost in a range between MinValue and MaxValue
make_kn_weighted(Size,MinValue,MaxValue,graph(LV,Comb)):-
    Size1 is Size+1,
    make_ordered_list(1,Size1,LV),
    find_all_combinations_weighted(LV,MinValue,MaxValue,[],Comb).

find_all_combinations_weighted([_],_,_,C,C):- !.
find_all_combinations_weighted([H|T],Min,Max,CT,CO):-
    find_combinations_weighted(H,T,Min,Max,C),
    append(CT,C,C1),
    find_all_combinations_weighted(T,Min,Max,C1,CO).

find_combinations_weighted(_,[],_,_,[]):- !.
find_combinations_weighted(E,[H|T],Min,Max,[edge(E,H,V)|TE]):-
    random(Min,Max,V),
    find_combinations_weighted(E,T,Min,Max,TE).


% make_kn_from_vertices(+Vertces,-Graph)
% creates a Kn Graph of size Size with given verteces name
make_kn_from_vertices(LV,graph(LV,Comb)):-
    find_all_combinations(LV,[],Comb).


% cycle_unweighted(+Graph,+Vert,-Cycle)
% find cycles in unweighted graph
cycle_unweighted(Graph,V,C):-
    find_path_unweighted(Graph,V,V,C),
    length(C,N),
    N > 3.


% cycle_weighted(+Graph,+Vert,-Cycle)
% find cycles in weighted graph
cycle_weighted(Graph,V,Cycle,Cost):-
    find_path_weighted(Graph,V,V,Cycle,Cost),
    length(Cycle,N),
    N > 3.    


% is_connected(+Graph)
% check if the graph is connected (checks for each vertex if there is a path
% to each remaining vertex)
% exist_path 0 for uneighted, 1 for weighted
is_connected(graph(LV,Edges)):-
    (memberchk(edge(_,_),(Edges)) -> 
        exist_path(LV,Edges,0);
    exist_path(LV,Edges,1)).

exist_path([_],_,_):-!.
exist_path([H|T],E,V):-
    exist_path_(H,T,E,V),
    exist_path(T,E,V).

exist_path_(_,[],_,_):-!.
exist_path_(H,[H1|T1],E,V):-
    (V = 0 ->
        connected(H,H1,E);
    connected_w(H,H1,E,_)),!,
    exist_path_(H,T1,E,V).


% node_degree(+Graph,?Node,?Deg).
% returns in Deg, the degree of Node
node_degree(graph(LN,Edges),Node,Deg):-
    member(Node,LN),
    (memberchk(edge(_,_),Edges) -> 
        (findall(X,(member(edge(X,Node),Edges)),VX),
            findall(Y,(member(edge(Node,Y),Edges)),VY)) 
        ;
        (findall(X,(member(edge(X,Node,_),Edges)),VX),
            findall(Y,(member(edge(Node,Y,_),Edges)),VY))
    ),
    length(VX,NX),
    length(VY,NY),
    Deg is NX+NY.


% node_degree_list(+Graph,-ListOfVerticesDegreeOrder)
% retuns a list of list of form [Degree,Edge]
node_degree_list(graph(LV,Edges),L):-
    node_degree_list_(LV,LV,Edges,[],L).

node_degree_list_([],_,_,L,Sorted):-
   sort(0,@>=,L,Sorted). % to keep equal values, msort
node_degree_list_([H|T],LV,Edges,LT,L):-
    node_degree(graph(LV,Edges),H,Deg),
    A = [Deg,H],
    append(LT,[A],LT1),
    node_degree_list_(T,LV,Edges,LT1,L).


% empty_unweighted_graph(+NumNodes,+NumEdges,-Graph)
% returns in Graph an empty graph 
empty_unweighted_graph(NumNodes,NumEdges,graph(LN,Edges)):-
    NumNodes > 0,
    NumEdges < NumNodes*(NumNodes - 1),
    length(LN,NumNodes),
    length(Edges,NumEdges),
    create_edges_unweighted(Edges).

create_edges_unweighted([]).
create_edges_unweighted([edge(_,_)|T]):-
    create_edges_unweighted(T).


% empty_weighted_graph(++NumNodes,++NumEdge,-Graph)
% returns in Graph an empty graph 
empty_weighted_graph(NumNodes,NumEdges,graph(LN,Edges)):-
    NumNodes > 0,    
    NumEdges < NumNodes*(NumNodes - 1),    
    length(LN,NumNodes),
    length(Edges,NumEdges),
    create_edges_weighted(Edges).

create_edges_weighted([]).
create_edges_weighted([edge(_,_,_)|T]):-
    create_edges_weighted(T).


% is_graph_node succeeds if a node is in the graph
% is_graph_node(+Graph,?Node).
% if Node is not ground, it returns all the nodes in the graph
is_graph_node(graph(LN,_),N):-
    (ground(N) ->
        memberchk(N,LN);
    member(N,LN)
    ).


% is_isolated_node succeeds if a node is isolated
% is_isolated_node(+Graph,?N)
% returns all isolated nodes in backtracking
is_isolated_node(graph(LN,Edges),N):-
    (memberchk(edge(_,_),Edges) ->
        (member(N,LN),
        \+member(edge(N,_),Edges),
        \+member(edge(_,N),Edges))
        ;
        (member(N,LN),
        \+member(edge(N,_,_),Edges),
        \+member(edge(_,N,_),Edges))
    ).


% is_graph_edge succeeds if a edge is in the graph
% is_graph_edge(+Graph,?Edge).
% if Edge is not ground, it returns all the nodes in the graph
is_graph_edge(graph(_,Edges),E):-
    (ground(E) ->
        memberchk(E,Edges);
    member(E,Edges)
    ).


% get_adjacent_nodes returns the adjacent nodes from the given one
% get_adjacent_nodes(+Graph,+Node,-List).
% if Node is not ground, returns all not nodes with corresponding list
get_adjacent_nodes(graph(LN,Edges),Node,List):-
    member(Node,LN),
    (memberchk(edge(_,_),Edges) ->
        (findall(X,(member(edge(Node,X),Edges)),VX),
        findall(Y,(member(edge(Y,Node),Edges)),VY))
        ;
        (findall(X,(member(edge(Node,X,_),Edges)),VX),
        findall(Y,(member(edge(Y,Node,_),Edges)),VY))
    ),
    append(VX,VY,VT),
    sort(VT,List).    


% reverse_edges(?Graph,?GraphRev)
% returns the graph with reversed edges
graph_reverse_edges(graph(LN,Edges),graph(LN,RevEdges)):-
    reverse_edges_(Edges,RevEdges).
reverse_edges_([],[]).
reverse_edges_([edge(X,Y,V)|T],[edge(Y,X,V)|T1]):-!,
    reverse_edges_(T,T1).
reverse_edges_([edge(X,Y)|T],[edge(Y,X)|T1]):-
    reverse_edges_(T,T1).


% spanning_trees(+Graph,-SpanningTree)
spanning_tree(graph([N|T],Edges),graph([N|T],TreeEdges)) :- 
   generate_spanning_tree(T,Edges,TreeEdgesUnsorted),
   sort(TreeEdgesUnsorted,TreeEdges).

generate_spanning_tree([],_,[]).
generate_spanning_tree(Curr,Edges,[Edge|T]) :- 
    select(Edge,Edges,Edges1), % select an edge and remove Edge from Edges
    get_vertices(Edge,X,Y), % find vertices adjacent
    is_connected_to_tree(X,Y,Curr), % check if connected
    delete(Curr,X,Curr1), % delete the two vertices
    delete(Curr1,Y,Curr2),
    generate_spanning_tree(Curr2,Edges1,T).

get_vertices(edge(X,Y),X,Y).
get_vertices(edge(X,Y,_),X,Y).

is_connected_to_tree(X,Y,Ns):- 
    memberchk(X,Ns), 
    \+ memberchk(Y,Ns), !.
is_connected_to_tree(X,Y,Ns):- 
    memberchk(Y,Ns), 
    \+ memberchk(X,Ns).


% mst_prim(+Graph,-MST)
% no choice poits left opened
mst_prim(graph([H|T],Edges),graph([H|T],TreeEdges),Cost):-
    predsort(compare_edges_value,Edges,SortedEdges),
    generate_spanning_tree(T,SortedEdges,TreeEdgesUnsorted),!, % keep it?
    sort(TreeEdgesUnsorted,TreeEdges),
    sum_cost(TreeEdges,0,Cost).

compare_edges_value(O,edge(X1,Y1,C1),edge(X2,Y2,C2)):-
    compare(O,C1+X1+Y1,C2+X2+Y2).

sum_cost([],C,C).
sum_cost([edge(_,_,C)|T],CT,Tot):-
    CT1 is CT+C,
    sum_cost(T,CT1,Tot).


% merge_graph(+Graph1,+Graph2,-MergedGraph)
% merges two graphs
% check if Edges1 and Edges2 have the same structure
merge_graph(graph(L1,Edges1),graph(L2,Edges2),graph(L,Edges)):-
    append(L1,L2,LT),
    sort(LT,L),
    append(Edges1,Edges2,EdgesT),
    sort(EdgesT,Edges).

% cut(+Graph,-G1,G2) generate a cut
% Finds a graph edge given source and target node
% finds all edges between source node and target node
% Returns a sorted list of all edges ending in TargetNode
% Returns the array of node names
% is acyclic
% graph modification: remove node, rename node, invert single edge
% flows algorithms Edmons Karp
% dual graph


% clique