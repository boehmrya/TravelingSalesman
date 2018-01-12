
/*
 * For code sources, I used this example for inspiration on kruskal:
 * http://www.geeksforgeeks.org/greedy-algorithms-set-2-kruskals-minimum-spanning-tree-mst/
 * I also used this example as inpiration for boruvka's:
 * http://www.geeksforgeeks.org/greedy-algorithms-set-9-boruvkas-algorithm/
 */

/**
 * @author Ryan Boehm
 */


import java.util.*;
import java.io.*;


public class TravelingSalesman {
    
    // Fields for Graph object
    public int V; // number of vertices in graph
    public int E; // number of edges in graph
    public LinkedList<Edge> edges = new LinkedList(); // all edges in graph
     
    // Static variable for cities array
    public static String[] cities;
    
    // Static variable for DFS vertices
    public static LinkedList DFSvertices = new LinkedList();
   
    
    // Creates a graph with V vertices and E edges
    hw8(int v, int e) {
        V = v;
        E = e;
        edges = new LinkedList();
    }
   
    
    // inner class to represent a graph edge
    class Edge implements Comparable<Edge> {
        
        // Fields for Edge object
        int source; // source vertex
        int dest; // destination vertex
        int weight; // distance in this case

        // Comparator function used for sorting edges based on
        // their weight
        public int compareTo(Edge otherEdge) {
            return this.weight - otherEdge.weight;
        }
        
        // Constructor for edges
        Edge(int u, int v, int w) {
            source = u;
            dest = v;
            weight = w; 
        }
        
        // default constructor
        Edge() {}
        
        // for printing out values of MST edges
        public String toString() {
            return "source: " + source + ", destination: " + dest + ", weight: " + weight;
        }
        
    }
    
    
    // setter method to add a new edge
    void addEdge(int u, int v, int w) {
        Edge thisEdge = new Edge(u, v, w);
        edges.add(thisEdge);
    }
    
 
    // inner class to represent a component of vertices/edges for union-find
    class component {
        int parent; // parent vertex of the component
        int rank; // used to compare smaller to larger components (depth)
        int size; // the number of nodes in the component
        
        // for printing out values of MST edges
        public String toString() {
            return " (parent: " + parent + ", rank: " + rank + ", size: " + size + ")";
        }
    };
   
 
    // A utility function to find component of a vertex i
    // (uses path compression technique)
    int find(component components[], int i) {
        // find root and make root as parent of i (path compression)
        if (components[i].parent != i) {
            components[i].parent = find(components, components[i].parent);
        }  
        return components[i].parent;
    }
    
    
    /*
    * calculates the smallest connected component for each component
    * that is still left.
    */
    public int smallestConnComp( component components[] ) {
        int minSize = 128; // set to max number possible initially
        for ( int i = 0; i < components.length; i++ ) {
            int thisSize = components[i].size;
            if ( thisSize > 0 ) {
                if ( thisSize < minSize ) {
                    minSize = thisSize;
                }
            }
        }
        return minSize;
    }
    
    
    /*
    * calculates the maximum depth of a reverse tree (see union method).
    */
    public int maxTreeDepth( component components[] ) {
        int maxDepth = 0;
        for ( int i = 0; i < components.length; i++ ) {
            int thisDepth = components[i].rank;
            if ( thisDepth > maxDepth ) {
                maxDepth = thisDepth;
            }
        }
        return maxDepth;  
    }
 
    
    // A function that does union of two compnents of x and y
    // (uses union by rank)
    // uses the shallow-threaded trees approach due to the find method
    // using path compression
    public void union(component components[], int x, int y) {
        int xroot = find(components, x);
        int yroot = find(components, y);
 
        // Attach smaller rank tree under root of high rank tree
        // (Union by Rank)
        if (components[xroot].rank < components[yroot].rank) {
            components[xroot].parent = yroot; // yroot parent of former xroot
            components[yroot].size = components[yroot].size + components[xroot].size;
            components[xroot].size = 0;
        }   
        else if (components[xroot].rank > components[yroot].rank) {
           components[yroot].parent = xroot;
           components[xroot].size = components[xroot].size + components[yroot].size;
           components[yroot].size = 0;
        }
 
        // If ranks are same, then make one as root and increment
        else {
            components[yroot].parent = xroot;
            components[xroot].size = components[xroot].size + components[yroot].size;
            components[yroot].size = 0;
            components[xroot].rank++; // stand-in for depth
        }
    }
    
    
    // The main function to construct MST using Boruvka's algorithm
    // returns the total weight of the MST
    // prints out each edge in the MST
    void BoruvkaMST() {            
 
        // Initially there are V different components (one for each vertex).
        // When the algorithm finishes, there will be one component that will be the MST
        int numComponents = V;
        int MSTweight = 0;
        Edge[] minWeightEdges = new Edge[V]; // hold min weight edges for each component
 
        // Create V components with single vertices
        component components[] = new component[V];
        for ( int i = 0; i < V; i++ ) {
            components[i] = new component();
            components[i].parent = i; // all parents are set to themsleves for one vertex components
            components[i].rank = 0; // all ranks set to 0 for one-vertex components
            components[i].size = 1; // set the size to 1 for all one-vertext components
            minWeightEdges[i] = new Edge(0,0,-1); // initialize all min weights to -1
        }
       
        // Keep combining components (or sets) until all
        // compnentes are not combined into single MST
        int numIterations = 0;
        
        while ( numComponents > 1 ) {
            
            // Traverse through all edges and update
            // minimum weight edge of every component
            for ( int i = 0; i < E; i++ ) {
 
                // Find components of two corners of current edge
                Edge thisEdge =  edges.get(i);
                int comp1 = this.find( components, thisEdge.source);
                int comp2 = this.find( components, thisEdge.dest);
 
                // If two vertices of current edge belong to
                // same component, ignore current edge. Else check if 
                // current edge is the minimum weight edge going out of the components
                // that it's endpoints u and v belong to.  Set the value in minWeightEdges if so.
                if ( comp1 != comp2 ) {
                    if ( (minWeightEdges[comp1].weight == -1) || ( minWeightEdges[comp1].weight > thisEdge.weight) ) {
                        minWeightEdges[comp1] = thisEdge;
                    } 
                    if ( (minWeightEdges[comp2].weight == -1) || ( minWeightEdges[comp2].weight > thisEdge.weight) ) {
                        minWeightEdges[comp2] = thisEdge;
                    } 
                }
                 
            }
            
            // Consider the above picked cheapest edges and add them to the MST
            for ( int i = 0; i < V; i++ ) {

                // Check if cheapest for current set exists
                if ( minWeightEdges[i].weight != -1 ) {
                    Edge safeEdge = minWeightEdges[i];
                    int comp1 = this.find( components, safeEdge.source );
                    int comp2 = this.find( components, safeEdge.dest );

                    if ( comp1 != comp2 ) {
                       MSTweight += safeEdge.weight;
                       this.union( components, comp1, comp2);
                       
                       // decrement number of components by one
                       numComponents--;                        
                    }
                }
            } 
            
            // increment the number of iterations to print out later
            numIterations++;
            
            // reset the min weight edges array
            Arrays.fill(minWeightEdges, new Edge(0,0,-1));
            System.out.println("Boruvka: Number of components: " + numComponents);
       
            // calculate the smallest connected component
            int smallestComp = smallestConnComp( components );
            System.out.println("Boruvka: Smallest Connected Component After " + numIterations + " iterations "
                    + "equals " + smallestComp);
            System.out.println();
        }
        
        // print out the total mst weight
        System.out.println("Boruvka: Number of Iterations: " + numIterations);
        System.out.println("Total Boruvka MST Weight: " + MSTweight);
        System.out.println();
        System.out.println();
    }
    
    
    // The main function to construct MST using Kruskal's algorithm
    Edge[] KruskalMST() {
        Edge result[] = new Edge[V];  //stores the resultant MST
        int MSTweight = 0;
        
        int e = 0;  // An index variable, used for result[]
        int i = 0;  // An index variable, used for sorted edges
        for ( i = 0; i < V; i++ ) {
            result[i] = new Edge();
        }
            
        // Step 1:  Sort all the edges in non-decreasing order of their
        // weight.  If we are not allowed to change the given graph, we
        // can create a copy of array of edges
        Collections.sort(edges);
 
        // Allocate memory for creating V ssubsets
        component components[] = new component[V];
        for( i = 0; i < V; i++ ) {
           components[i] = new component(); 
        }
            
        // Create V subsets with single elements
        for ( int j = 0; j < V; j++ ) {
            components[j].parent = j;
            components[j].rank = 0;
        }
 
        i = 0;  // Index used to pick next edge
 
        // Number of edges to be taken is equal to V-1
        while (e < V - 1) {
            // Step 2: Pick the smallest edge. 
            // increment the index for next iteration
            Edge nextEdge = edges.get(i);
            i++;
 
            int comp1 = find(components, nextEdge.source);
            int comp2 = find(components, nextEdge.dest);
 
            // If including this edge does't cause cycle, include it
            // in result and increment the index of result for next edge
            if (comp1 != comp2) {
                result[e++] = nextEdge;
                MSTweight += nextEdge.weight;
                union(components, comp1, comp2);
                
                //calculate the depth of the maximum reverse tree.
                int maxDepth = maxTreeDepth( components );
                System.out.println("maximum tree depth after iteration " + e + " is " + maxDepth);
            }
            // otherwise just pass over the edge (discard) 
        }
 
        // print the contents of result[] to display the built MST
        System.out.println(); // line break
        System.out.println("Following are the edges in the Kruskal MST");
        for (i = 0; i < e; ++i) {
            String srcName = cities[result[i].source];
            String dstName = cities[result[i].dest];
            System.out.println( srcName + " to " + dstName + " with weight " +
                               result[i].weight);
        }
        System.out.println(); // line break
        System.out.println("Kruskal MST Weight = " + MSTweight);
        System.out.println(); // line break
        return result;
    }
 

    /*
     * method to count the number of cities in the file to initialize the adjacency matrix.
     */
    public static int countCities( String filename ) {
        String line;
        int numCities = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader( filename ));
            
            // go line by line and count the number of cities
            while (((line = br.readLine()) != null)){
                char c = line.charAt(0);
                if ((c != '*') && (Character.isLetter(c))){
                    numCities++;         
                }
            }
            
            // close the file
            br.close();
            
        } catch(IOException e){
            System.out.println( "File not found" );
        }
        
        return numCities;
    }
    
    
    /*
     * method to count the number of edges in the file (graph).
     */
    public static int countEdges( String filename ) {
        String line;
        int numEdges = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader( filename ));
            
            // go line by line and count the number of edges
            while (((line = br.readLine()) != null)){
                char c = line.charAt(0);
                if ((c != '*') && (Character.isDigit(c))) {
                    String[] thisLineEdges = line.split("\\s+");
                    for ( int i = 0; i < thisLineEdges.length; i++ ) {
                        numEdges += 1;
                    }
                }
            }
            
            // close the file
            br.close();
            
        } catch(IOException e){
            System.out.println( "File not found" );
        }
        
        return numEdges;
    }
    
    
    /*
    * Build the cities array with the name of each city
    */
    public static String[] InitCities( String fileName, int numCities ){
        String line;
        int cityIndex = 0;
        String[] cities = new String[numCities];
        
        try {
            BufferedReader br = new BufferedReader( new FileReader( fileName ));
            
            while (((line = br.readLine()) != null)){
                char c = line.charAt(0);
                if ((c != '*') && (Character.isLetter(c))){
                    if ((c != '*')){
                        String[] cityInfo = line.split(",");
                        String thisCity = cityInfo[0];
                        cities[cityIndex] = thisCity;
                        cityIndex++;
                    } 
                }  
            }
            
            br.close();  
        } catch(IOException e){
            System.out.println(e);
        }
        return cities;
    }
    
    
    /*
    * Builds the graph object by initializing the object and creating all
    * corresponding edges.
    */
    public static hw8 InitGraph( String fileName, int numCities, int numEdges ){
        String line;
        hw8 thisGraph = new hw8(numCities, numEdges);
        int cityCounter = -1;
        int cityEdgeCounter = 1;
    
        try {
            BufferedReader br = new BufferedReader( new FileReader( fileName ));
            
            while (((line = br.readLine()) != null)){
                char c = line.charAt(0);
                // line is a city
                if ((c != '*') && (Character.isLetter(c))){
                    cityCounter++;
                    cityEdgeCounter = 1; // reset city edge counter to start for next city
                } 
                // line has edges
                else if ((c != '*') && (Character.isDigit(c))) {
                    String[] thisLineEdges = line.split("\\s+");
                    // go through each edge, create edge object, add to graph
                    for ( int i = 0; i < thisLineEdges.length; i++ ) {
                        int source = cityCounter;
                        int dest = cityCounter - cityEdgeCounter;
                        int weight = Integer.parseInt(thisLineEdges[i]);
                        thisGraph.addEdge(source, dest, weight);
                        cityEdgeCounter++;
                    }
                    
                }
            }
            
            br.close();  
        } catch(IOException e){
            System.out.println(e);
        }
        return thisGraph;
    }
    
    
    /*
    * Takes an array of edge objects and converts
    * it into an adjacency list, represented by an 
    * array of linked lists.
    */
    LinkedList[] createAdjList( Edge edges[] ) {

        // intilize the adjacency list
        LinkedList[] adjList = new LinkedList[V];
        for (int i = 0; i < V; i++ ) {
            adjList[i] = new LinkedList();
        }
        
        // add all vertices to the adjacency list
        for( int i = 0; i < edges.length; i++ ) {
            int u = edges[i].source;
            int v = edges[i].dest; 
            
            // add vertices to adjacency list
            if (!adjList[u].contains(v)) {
               adjList[u].add(v); 
            }
            if (!adjList[v].contains(u)) {
                adjList[v].add(u);
            } 
        }   
        return adjList;
    }
    
    
    
    /*
    * Utility function used by Depth-First Search to do recursive calls
    */
    void DFSUtil(int v, boolean visited[], LinkedList[] adjList) {
        // Mark the current node as visited and print it
        visited[v] = true;
        DFSvertices.add(v);
 
        // Recur for all the vertices adjacent to this vertex
        Iterator<Integer> i = adjList[v].listIterator();
        while (i.hasNext()) {
            int n = i.next();
            if (!visited[n]) {
               DFSUtil(n, visited, adjList); 
            }    
        }
    }
 
    
    /* 
    * The function to do DFS traversal. It uses recursive DFSUtil()
    */
    void DFS(int v, LinkedList[] adjList) {
        // Mark all the vertices as not visited(set as
        // false by default in java)
        boolean visited[] = new boolean[this.V];
 
        // Call the recursive helper function to print DFS traversal
        DFSUtil(v, visited, adjList);
    }
    
    
    /*
    * takes a linked list of vertices and removes duplicates
    */
    LinkedList<Integer> removeDupes(LinkedList<Integer> vertices) {
        LinkedList finalVertices = new LinkedList();
        for (int i = 0; i < vertices.size(); i++ ) {
            int thisVertex = vertices.get(i);
            if (!finalVertices.contains(thisVertex)) {
                finalVertices.add(thisVertex);
            }
        }
        return finalVertices;   
    }
    
    
    /*
    * takes as input an array of edges for an MST
    * and calculates the metric TSP and returns an array of edge objects
    */
    LinkedList<Edge> metricTSP( Edge[] mstEdges ) {
        // create adjacency list of MST result (for DFS below)
        LinkedList[] mstAdjList = this.createAdjList( mstEdges );
        LinkedList<Edge> TSPedges = new LinkedList();
                
        // depth first search starting at vertext
        // adds raw duplicated vertices to DFSvertices list field in graph class
        this.DFS(0, mstAdjList);
        
        // remove dupes from vertices after DFS
        LinkedList<Integer> TSPvertices = this.removeDupes( DFSvertices );
                
        // go through each pair of vertices in the list
        for (int i = 0; i < TSPvertices.size() - 1; i++ ) {
            int u = TSPvertices.get(i);
            int v = TSPvertices.get(i + 1);
            
            
            // get the corresponding edge (source, dest, and weight) for 
            // all of the edges in the TSP except the last
            for (int j = 0; j < edges.size() - 1; j++ ) {
                Edge thisEdge = edges.get(j);
                
                // determine if either source or destination of edge match
                // the vertices in the TSP
                if ((thisEdge.source == u) && (thisEdge.dest == v)) {
                    Edge thisTSPedge = new Edge(u, v, thisEdge.weight);
                    TSPedges.add(thisTSPedge);     
                }
                else if ((thisEdge.source == v) && (thisEdge.dest == u)) {
                    Edge thisTSPedge = new Edge(u, v, thisEdge.weight);
                    TSPedges.add(thisTSPedge); 
                }
            }  
        }
        
        // add the final leg of the TSP
        int u = TSPvertices.get(TSPvertices.size() - 1); // the last vertex
        int v = TSPvertices.get(0); // the beginning vertex
        
        // find the edge that corresponds to the last leg in the TSP
        for (int j = 0; j < edges.size() - 1; j++ ) {
            Edge thisEdge = edges.get(j);
                
            // determine if either source or destination of edge match
            // the vertices in the TSP
            if ((thisEdge.source == u) && (thisEdge.dest == v)) {
                Edge thisTSPedge = new Edge(u, v, thisEdge.weight);
                TSPedges.add(thisTSPedge);     
            }
            else if ((thisEdge.source == v) && (thisEdge.dest == u)) {
                Edge thisTSPedge = new Edge(u, v, thisEdge.weight);
                TSPedges.add(thisTSPedge); 
            }
        }
        
        return TSPedges;
    }
    
    
    /*
    * Takes list of all edges in the TSP and returns the total weight;
    */
    int metricTSPWeight( LinkedList<Edge> TSPedges) {
        int TSPweight = 0;
        for (int i = 0; i < TSPedges.size(); i++ ) {
            Edge thisEdge = TSPedges.get(i);         
            int thisEdgeWeight = thisEdge.weight;
            TSPweight += thisEdgeWeight;
        }
        return TSPweight;
    }
    
    
    /*
    * check if the tsp is a valid tsp - used after swaps
    * to check whether swap should be retained.
    */
    boolean isValidTSP ( LinkedList<Edge> TSPedges ) {
        // check the number of edges first
        if (TSPedges.size() != 128) {
            return false;
        }
       // determine if the source of the next edge equals the
       // destination of the previous edge.  If this is not the case
       // in a single instance, return false.  Otherwise, return true.
       for (int i = 0; i < (TSPedges.size() - 1); i++ ) {
           int u = TSPedges.get(i).source;
           int v = TSPedges.get(i).dest; 
           
           // find any duplicate source or destination edges
           for (int j  = i + 1; j < TSPedges.size(); j++ ) {
               int newU = TSPedges.get(j).source;
               int newV = TSPedges.get(j).dest;
               
               // if sources or destinations are duplicated after swap
               // then return false
               if ( (u == newU) || (v == newV) ) {
                   return false;
               }
               
           } 
       }
       return true;
    }
    
    
    
    /*
    * Two-opt traveling salesman tour
    */
    LinkedList<Edge> opt2TSP( LinkedList<Edge> TSPedges, int TSPweight ) {
        
        // the new TSP to return
        LinkedList<Edge> newTSP = TSPedges;
        int newTSPweight = TSPweight;
        int edgeSwap = 0; // to track whether swaps occurred 
        
        // for each edge in the tsp
        for (int i = 0; i < (TSPedges.size() - 2); i++ ) {
            
            // first edge u and v values
            Edge firstEdge = TSPedges.get(i);
            int firstU = firstEdge.source;
            int firstV = firstEdge.dest;
            
            // check every other remaining edge for potential swaps
            for (int j = i + 2; j < TSPedges.size(); j++ ) {
                
                // second edge u and v values
                Edge secondEdge = TSPedges.get(j);
                int secondU = secondEdge.source;
                int secondV = secondEdge.dest;
                
                // to track actual edge swaps
                // set to 0 before going through all edges.
                edgeSwap = 0;
                
                // copy the linked list of the TSP
                LinkedList tempTSP = new LinkedList();
                tempTSP = (LinkedList) newTSP.clone();
                
                // find two edges you can swap out
                for (int k = 0; k < edges.size(); k++ ) {
                    Edge thisEdge = edges.get(k);
                    // ensure that the tour doesn't contain the reverse of this edge
                    Edge reverseEdge = new Edge(thisEdge.dest, thisEdge.source, thisEdge.weight);
                    
                    // only swap edges if this edge is not the reverse of a current edge in the TSP
                    if (tempTSP.contains(reverseEdge) == false) {
                        if ((thisEdge.source == firstU) && (thisEdge.dest == secondV)) {
                            tempTSP.set(i, thisEdge);
                            edgeSwap++;
                        }
                        else if ((thisEdge.source == secondU) && (thisEdge.dest == firstV)) {
                            tempTSP.set(j, thisEdge);
                            edgeSwap++;
                        }
                    }
                    
                    
                    // if you've made two swaps, break out of loop
                    if ( edgeSwap == 2 ) {
                        
                        // validity and weight
                        boolean isValid = isValidTSP(tempTSP);
                        int tempTSPweight = metricTSPWeight(tempTSP);
                        
                        if ( (isValid == true) && (tempTSPweight < newTSPweight) )  {
                            newTSP = tempTSP;
                        }

                        // reset edge swap counter and break
                        edgeSwap = 0;
                        break;
                    }
                }
                  
            }  
        }
        return newTSP;
    }
    
   
    
    
    /*
    * Three-opt traveling salesman tour
    */
    LinkedList<Edge> opt3TSP( LinkedList<Edge> TSPedges, int TSPweight ) {
        // the new TSP to return
        LinkedList<Edge> newTSP = TSPedges;
        int newTSPweight = TSPweight;
        int edgeSwap = 0; // to track whether swaps occurred 
        
        // for each edge in the tsp
        for (int i = 0; i < (TSPedges.size() - 4); i++ ) {
            
            // first edge u and v values
            Edge firstEdge = TSPedges.get(i);
            int firstU = firstEdge.source;
            int firstV = firstEdge.dest;
            
            // get a second edge for potentials swaps
            for (int j = i + 2; j < TSPedges.size() - 2; j++ ) {
                
                // second edge u and v values
                Edge secondEdge = TSPedges.get(j);
                int secondU = secondEdge.source;
                int secondV = secondEdge.dest;
                
                // get a third edge for potential swaps
                for (int h = j + 2; h < TSPedges.size(); h++ ) {
                    
                    // second edge u and v values
                    Edge thirdEdge = TSPedges.get(j);
                    int thirdU = thirdEdge.source;
                    int thirdV = thirdEdge.dest;
                    
                    // to track actual edge swaps
                    // set to 0 before going through all edges.
                    edgeSwap = 0;

                    // copy the linked list of the TSP
                    LinkedList tempTSP = new LinkedList();
                    tempTSP = (LinkedList) newTSP.clone();

                    // find two edges you can swap out
                    for (int k = 0; k < edges.size(); k++ ) {
                        Edge thisEdge = edges.get(k);
                        if ((thisEdge.source == firstU) && (thisEdge.dest == secondV)) {
                            tempTSP.set(i, thisEdge);
                            edgeSwap++;
                        }
                        if ((thisEdge.source == firstU) && (thisEdge.dest == thirdV)) {
                            tempTSP.set(i, thisEdge);
                            edgeSwap++;
                        }
                        else if ((thisEdge.source == secondU) && (thisEdge.dest == firstV)) {
                            tempTSP.set(j, thisEdge);
                            edgeSwap++;
                        }
                        else if ((thisEdge.source == secondU) && (thisEdge.dest == thirdV)) {
                            tempTSP.set(j, thisEdge);
                            edgeSwap++;
                        }
                        else if ((thisEdge.source == thirdU) && (thisEdge.dest == firstV)) {
                            tempTSP.set(h, thisEdge);
                            edgeSwap++;
                        }
                        else if ((thisEdge.source == thirdU) && (thisEdge.dest == secondV)) {
                            tempTSP.set(h, thisEdge);
                            edgeSwap++;
                        }

                        // if you've made two swaps, break out of loop
                        if ( edgeSwap == 3 ) {

                            // validity and weight
                            boolean isValid = isValidTSP(tempTSP);
                            int tempTSPweight = metricTSPWeight(tempTSP);

                            if ( (isValid == true) && (tempTSPweight < newTSPweight) )  {
                                newTSP = tempTSP;
                            }

                            // reset edge swap counter and break
                            edgeSwap = 0;
                            break;
                        }
                    }
                    
                }
      
            }  
        }
        return newTSP;  
    }
    
    /*
    * our main program
    */
    public static void main(String[] args) {
        String fileName = "miles.dat";
        
        // calculate the number of vertices (cities) and edges
        int numCities = countCities( fileName );
        int numEdges = countEdges( fileName );
                
        // fill out the cities array: a length 128 array with each
        // index holding a name of a city
        cities = InitCities( fileName, numCities);
        
        // build the graph with all edge objects
        hw8 roadGraph = InitGraph( fileName, numCities, numEdges);
        
        // run boruvka's MST algorithm on the graph
        //componentsLeft = new int[numCities];
        //Arrays.fill(componentsLeft, 1); // initialize all components to 1 (set to 0 as they are merged)
        System.out.println("Boruvka's MST: ");
        roadGraph.BoruvkaMST();
        
        // two spaces to separate the MST algorithms
        System.out.println();
        System.out.println();
        
        // run kruskal's MST algorithm on the graph
        System.out.println("Krustkal's MST: ");
        Edge[] mstEdges = roadGraph.KruskalMST();
                
        // two spaces to separate the MST algorithms
        System.out.println();
        System.out.println();
        
        // calculate the metric tsp using kruskal's mst
        LinkedList<Edge> metricTSP = roadGraph.metricTSP(mstEdges);
        int TSPWeight = roadGraph.metricTSPWeight(metricTSP);
        System.out.println("Total Metric TSP Weight: " + TSPWeight);
        
        // two spaces to separate the TSP algorithms
        System.out.println();
        System.out.println();
        
        // print the contents of result[] to display the built MST
        System.out.println("Following are the edges in the metric TSP");
        for (int i = 0; i < metricTSP.size(); i++) {
            String srcName = cities[metricTSP.get(i).source];
            String dstName = cities[metricTSP.get(i).dest];
            System.out.println( srcName + " to " + dstName + " with weight " +
                               metricTSP.get(i).weight);
        }
        
        // two spaces to separate the TSP algorithms
        System.out.println();
        System.out.println();
        
        // calculate the two-opt tsp
        LinkedList<Edge> opt2tsp = roadGraph.opt2TSP( metricTSP, TSPWeight);
        int opt2TSPWeight = roadGraph.metricTSPWeight(opt2tsp);
        System.out.println("Total Opt-2 TSP Weight: " + opt2TSPWeight);
         
        // line breaks
        System.out.println();
        System.out.println();
        
                 
         // print out the edges in the 2-opt solution
        System.out.println("Following are the edges in the Opt-2 TSP");
        for (int i = 0; i < opt2tsp.size(); i++) {
            String srcName = cities[opt2tsp.get(i).source];
            String dstName = cities[opt2tsp.get(i).dest];
            System.out.println( srcName + " to " + dstName + " with weight " +
                               opt2tsp.get(i).weight);
        }
       
        // two spaces to separate the TSP algorithms
        System.out.println();
        System.out.println();
        
        // calculate the three-opt tsp
        LinkedList<Edge> opt3tsp = roadGraph.opt2TSP( opt2tsp, opt2TSPWeight);
        int opt3TSPWeight = roadGraph.metricTSPWeight(opt3tsp);
        System.out.println("Total Opt-3 TSP Weight: " + opt3TSPWeight);
        
        // two spaces to separate the TSP algorithms
        System.out.println();
        System.out.println();
        
        // print out the edges in the 3-opt solution
        System.out.println("Following are the edges in the Opt-3 TSP");
        for (int i = 0; i < opt3tsp.size(); i++) {
            String srcName = cities[opt3tsp.get(i).source];
            String dstName = cities[opt3tsp.get(i).dest];
            System.out.println( srcName + " to " + dstName + " with weight " +
                               opt3tsp.get(i).weight);
        }
        
        
    }
    

}
