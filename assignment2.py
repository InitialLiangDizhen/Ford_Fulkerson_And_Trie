"""
FIT2004 Assignment 2 - 2023
Student Name: LIANG DIZHEN
Student ID: 31240291
Version 3.2.2
"""

from collections import deque
class Vertex:
    """
    Writtern by LIANG DIZHEN
    Vertex class for flow network and trie are modified from the vertex class in my assignment 1
    """

    def __init__(self, id):
        """
        The vertex object is created specialise for flow network and trie function which has vertex id, its outgoing edge
        status of discovered and visited,
        previous: edge that from another vertex to this vertex

        Input:
        id: vertex id


        Return:
            No return

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(1)
        """

        # vertex id
        self.id = id
        # edges list
        self.edges = []

        # for bfs
        self.discovered = False
        self.visited = False

        # backtracking/where I was from
        self.previous = None

    def __str__(self):
        """
        Write by LIANG DIZHEN
        this function is modified from the str function in vertex class in my assignment 1
        function to return basic information about this vertex object

        Input: No Input

        Return: string that has the basic infomation about the vertex

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(E) where E is the number of edge inside this vertex object
        """
        # print self.id
        return_string = str(self.id)

        # append its edge into a return string to print out
        for i in range(len(self.edges)):
            return_string = return_string + "\n with edges: " + str(self.edges[i])
        return return_string


    def add_edge(self, edge):
        """
        Write by LIANG DIZHEN
        function to add edge to the vertex object
        Input: edge object

        Return: No Return

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(1)
        """
        self.edges.append(edge)


class Edge:
    def __init__(self, u, v, c):
        """
        Write by LIANG DIZHEN
        Construct Edge object for flow network and trie

        Input: u - start vertex id, v - end vertex id, c - capacity of the edge

        Return: No Return

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(1)
        """

        self.u = u
        self.v = v
        self.c = c

        #re flow of edge in residual graph
        self.flow = 0
        self.previous = None

    def __str__(self):
        """
        Written by LIANG DIZHEN
        function to return basic information about this edge object

        Input: No input

        Return: string of basic information about this edge object

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(1)
        """
        return_string = str(self.u) + "," + str(self.v) + "," + " capacity: " + str(self.c) + " flow: " + str(self.flow)
        return return_string

    #check if edge has remaining capacity
    def hasCapacity(self):
        """
        Written by LIANG DIZHEN
        function to check if edge has remaining capacity

        Input: No input

        Return: true if edge has remaining capacity, false if edge has no remaining capacity

        Time Complexity:
        Best: O(1)
        Worst: O(1)

        Space Complexity:
        Input: O(1)
        Aux: O(1)
        """

        if (self.c - self.flow) > 0:
            return True
        return False

class Graph:
    """
    Writtern by LIANG DIZHEN
    This class is modified from the graph class in my assignment 1
    class for both question 1 and question2 to construct flow network and trie
    since both of them need to construct a graph
    """


    def __init__(self, argv_vertices_count):
        """
        Written by LIANG DIZHEN
        function to construct Graph object for both question 1 and question2 with the aid of add_Connection function
        Itself just construct a list of vertices for the graph

        Input: Number of vertices in the graph

        Return: No return

        Time Complexity:
        Best: O(V)
        Worst: O(V) V is the number of vertices in the graph

        Space Complexity:
        Input: O(1)
        Aux: O(V) where V is the number of vertices in graph
        """

        # start index for intermediate_1 vertices
        self.inter_1_start = argv_vertices_count + 1

        #start index for intermediate_2 vertices
        self.inter_2_start = self.inter_1_start * 2

        # total number of vertices - actual vertices + residual vertices
        #index = len(argv_vertices_count) max_vertices[index] = imaginary end
        #max_vertices =  argv_vertices_count + 2 = original vertices + imaginary end + imaginary start
        # *2 - 1 = intermediate vertex in between all vertices
        self.max_vertices = (argv_vertices_count * 3) + 2

        # imaginary start index
        self.imag_start_index = argv_vertices_count

        # imaginary end index
        self.imag_end_index = argv_vertices_count + self.inter_1_start


        # Array to construct Graph
        self.vertices = [None] * self.max_vertices  # list of vertices
        # create vertex object to store in Graph array
        for i in range(self.max_vertices):
            # create vertex object and add into the vertices list
            v = Vertex(i)
            self.vertices[i] = v

    def __str__(self):
        """
        Written by LIANG DIZHEN
        function to print basic information about the graph

        Input: No input

        Return: string contain basic information about the graph

        Time Complexity:
        Best: O(1)
        Worst: O(V) V is the number of vertex objects in the graph

        Space Complexity:
        Input: O(1)
        Aux: O(V)
        """

        return_string = ""
        for vertex in self.vertices:
            return_string = return_string + "Vertex" + str(vertex) + "\n"
        return return_string


    def add_Connection(self, connections, maxIn, maxOut, targets):
        """
        Writtern by LIANG DIZHEN
        This class is modfiied from the add_Road function in my assignment 1

        function to create edge and add into the corresponding vertex from the connections in the FlowNetwork
        create connections between intermediate_1 vertices and original vertices,
        create connections between original vertices and intermediate_2 vertices,
        create connections between u intermediate_2 vertices and v intermediate_1 vertices

        intermediate_1 vertices are used to for maxIn to limit the flow from incoming edges
        intermediate_2 vertices are used to for maxOut to limit the flow from outgoing edges

        Input: connections - all the possible connections between original nodes

        Return: no return

        Time Complexity:
        Best: O|D| + |C|)
        Worst: O(|D| + |C|)  D is the number of data centres, C is the number of connections

        Space Complexity:
        Input: O(|D| + |C|)
        Aux: O(|D| + |C|)
        """
        #create all intermediate imaginary vertices two for given each given vertex O(|D| + |C|)
        #one is for maxIn, one is for maxOut
        #Edge(v_id_1, v_id_2, capacity)
        #create for maxIn
        for i in range(len(maxIn)):
            #create intermediate vertex
            #maxIn from u intermediate_1 vertex to u vertex
            inter1_vid = i + self.inter_1_start
            inter1_vertex = self.vertices[inter1_vid]
            #create edge from intermediate 1 vertex to i vertex
            current_edge = Edge(inter1_vid, i, maxIn[i])
            inter1_vertex.add_edge(current_edge)
            #print("inter1_vid: " + str(inter1_vid) + " i: " + str(i) + " maxIn[i]: " + str(maxIn[i]))

        #create edge for maxOut
        #maxOut from u vertex to u_intermediate_2 vertex
        for i in range(len(maxOut)):
            #create intermediate2 vertex
            inter2_vid = i + self.inter_2_start
            #create edge from i vertex to intermediate2 vertex
            current_edge = Edge(i, inter2_vid, maxOut[i])
            #for the targets (actual end) only cares about maxInt
            if i in targets:
                current_edge = Edge(i, inter2_vid, maxIn[i])
            self.vertices[i].add_edge(current_edge)
            #print(" i: " + str(i) + " inter2_vid: " + str(inter2_vid) + " maxOut[i]: " + str(maxOut[i]))

        #single source, accumulate all outgoing capacity for start vertex
        start_num = 0

        #accmulate corresponding capacity of incoming edges to end vertex
        end_num = len(targets)
        # accumulative values in list for end vertex
        acc_lst = [0] * end_num

        #create edge from u inter_2 vertex to v inter_1 vertex O(|C|)
        for road in connections:
            u, v, c = road

            # add outgoing edge from u intermediate2 to v intermediate 1
            #intermediate 2 of u
            int_2_uid = u + self.inter_2_start
            int_2_u = self.vertices[int_2_uid]

            #intermediate 1 of v
            int_1_vid = v + self.inter_1_start

            #add Edge in between
            normal_edge = Edge(int_2_uid, int_1_vid, c)
            int_2_u.add_edge(normal_edge)

            #accumulative outgoing capacity to start vertex
            # if u == origin:
            #     start_num += c
                #print("start num: " + str(start_num))

            # accumulative incoming capacity to each end vertex
            # for i in range(end_num):
            #     if v == targets[i]:
            #         acc_lst[i] += c

        #create edge from imaginary start to intermediate1 vertex of u
        # img_start_v = self.vertices[self.imag_start_index]
        # origin_inter1_vid = origin + self.inter_1_start
        # img_start_edge = Edge(self.imag_start_index, origin_inter1_vid, maxIn[origin])
        # img_start_v.add_edge(img_start_edge)

        #create edge from intermediate2 vertex for targets of v to imaginary end
        for i in range(end_num): #O(|D|)
            # the vertices in targets is the actual end, no need to care about the
            # maxOut, even when it now connect to imaginary end
            # the targets vertices are just limited by the maxIn
            #print("targets: " + str(targets[i]))
            ver_out = targets[i] + self.inter_2_start
            end_edge = Edge(ver_out, self.imag_end_index, maxIn[targets[i]])


            current_vertex = self.vertices[ver_out]
            current_vertex.add_edge(end_edge)



class FlowNetwork(Graph):
    def __init__(self, argv_vertices_count):
        """
        Writtern by LIANG DIZHEN

        This init function use the super class Graph's init function to create a flow network
        with the given number of vertices, read graph init function for more information

        Precondition:
        Postcondition:

        Input:
            argv_vertices_count: number of data centres
        Return:
            flow network

        Time Complexity:
        Best: O(D) D is the number of data centres
        Worst: O(D) D is the number of data centres

        Space Complexity:
        Input: O(D)
        Aux: O(D)
        """
        fn = super().__init__(argv_vertices_count)

        return fn

    def reset(self):
        """
        Writtern by LIANG DIZHEN

        this function is modified from the super class Graph's reset function in my assingment1
        this function reset all vertices in the graph to undiscovered, unvisited and previous to None

        Input:
            None
        Return:
            None

        Time complexity:
            Best: O(V) where V is the number of vertices in the graph
            Worst: O(V) where V is the number of vertices in the graph
        Space complexity:
            Input: O(1)
            Aux: O(1)
        """
        for vertex in self.vertices:
            vertex.discovered = False
            vertex.visited = False
            vertex.previous = None


    def bfs(self, source, target):
        """
        Writtern by LIANG DIZHEN

        bfs is a function that use breadth first search to find a path from source to target
        in the flow network to find a path with overall residual capacity greater than 0
        and set the previous of each vertex to the previous edge that connected from another vertex to this vertex
        in the path

        Input:
            source: the source vertex id
            target: the target vertices id
        Return:
            True if there is a path from source to target
            False if there is no path from source to target

        Time complexity:
            Best: O(|C| * |D|) where C is the number of edges in the graph
            Worst: O(|C| * |D|) where D is the number of data centres
        Space complexity:
            Input: O(|D|) where D is the number of data centres
            Aux: O(|C| + |D|)
        """
        #since running multiple time
        self.reset()  # reset all vertices to undiscovered and unvisited
        source = self.vertices[source]
        discovered = deque()
        discovered.append(source)

        #while there is still vertex in the deque
        while len(discovered) > 0:
            u = discovered.popleft()  # serve() with complexity O(1)
            #print(u.id)
            u.discovered = True
            u.visited = True

            if u.id == target:
                return True

           #print("u_id: " + str(u))
            for edge in u.edges:
                v = edge.v  # v is the vertex id
                v = self.vertices[v]

                if edge.hasCapacity():
                    if v.discovered == False and v.visited == False:
                        discovered.append(v)
                        v.discovered = True
                        #get the edge then travel it to v
                        v.previous = edge

        return False # no path found


    def ford_fulkerson(self, s, t):
        """
        Written by LIANG DIZHEN

        ford_fulkerson is a function that use the bfs to find  max possible flow
        in the augmenting path that is discovered by the bfs, and the same time which would set the edge
        to be the previous of the destination vertex. Then use the max residual capacity augment the flow
        alone the path during backtracking by using previous of each vertex. After that, add a reverse edge from destination
        vertex to start vertex with the max possible flow to be its capacity.Then repeat the process until there is no
        augmenting path found by bfs.


        Precondition:
        Postcondition:

        Input:
            s - origin vertex id, t - target vertex id
        Return:
            max_flow - the accumulative flow of the max possible flow of each augmenting path found by bfs
            max flow that is possible to flow from s to t

        Time complexity:
            Best: O(|D| * |C|^2) where C is the number of edges in the graph
            Worst: O(|D| * |C|^2) where D is the number of data centres

        Space complexity:
            Input: O(|D|) where D is the number of data centres
            Aux: O(|C| + |D|)
        """

        #initialize max_flow to 0
        max_flow = 0

        #print("go int bfs")

        #while there is a path from s to t in the residual graph (find Augmenting Path)

        while self.bfs(s, t):
            #there is a augmenting path when get into this loop
            #v.previou must != None
            #to get the lowest flow in the residual graph
            path_flow = float("Inf")

            v = self.vertices[t]

            #stop when u.id = s
            while v.id != s:
                current_edge = v.previous
                u_id = current_edge.u
                u = self.vertices[u_id]

                possible_flow = current_edge.c - current_edge.flow

                if possible_flow < path_flow:
                    path_flow = possible_flow

                #update v until s
                v = u

            #run again to augment flow
            #update residual graph (remaining capacity by augmenting used flow on edge.flow)
            v = self.vertices[t]
            while v.id != s:
                current_edge = v.previous
                u = current_edge.u
                u = self.vertices[u]

                #augument the flow on forward edge
                current_edge.flow = current_edge.flow + path_flow

                    #adjuct back to true vertex id
                    # print("u: " + str(current_edge.u - self.inter_2_start)
                    #       + " v: " + str(current_edge.v - self.inter_1_start)
                    #       + " capacity: " + str(current_edge.c)
                    #       + " flow: " + str(current_edge.flow))
                    # current_edge.flow = 0

                # #add reverse edge from v intermediate 2 to u intermediate 1
                rev_edge = Edge(v.id, u.id, path_flow)
                rev_edge.flow = 0

                #add reverse edge from v to u to allow backward flow
                #or flow to be undone
                v.add_edge(rev_edge)
                # backtrack to u
                v = u

            #print("path_flow", max_flow)
            max_flow += path_flow

        return max_flow


def maxThroughput(connections, maxIn, maxOut, origin, targets):
    """
    Written by LIANG DIZHEN
    maxThroughput would build the flow network first by creating a list to store all vertices
    (original vertices + intermediate 1 vertex to limit the flow from incoming edge for each vertex,
    intermediate 2 vertex to limit the flow of outgoing edge for each vertex)

    After that, Build the flow network by adding edges to the corresponding vertices in the list.
    Add edge from intermediate 1 vertex to original vertex with capacity of maxIn for each vertex
    Add edge from original vertex to intermediate 2 vertex with capacity of maxOut for each vertex
    Add edge from u intermediate 2 vertex to v intermediate 1 vertex with capacity of connection from u to v

    After that, use ford_fulkerson to find the max possible flow by going through all the augmenting path with bfs

    ford_fulkerson is a function that use the bfs to find  max possible flow
    in the augmenting path that is discovered by the bfs, and the same time which would set the edge
    to be the previous of the destination vertex. Then use the max residual capacity augment the flow
    alone the path during backtracking by using previous of each vertex. After that, add a reverse edge from destination
    vertex to start vertex with the max possible flow to be its capacity.Then repeat the process until there is no
    augmenting path found by bfs.


    Precondition:
    Postcondition:

    Input:
        maxIn - list of max possible flow into each data centre
        maxOut - list of max possible flow out of each data centre
        connections - list of connections between data centres
        origin - the origin vertex id
        targets - list of target vertex id

    Return:
        max_flow - the accumulative flow of the max possible flow of each augmenting path found by bfs
        max flow that is possible to flow from s to t

    Time complexity:
        Best: O(|D| * |C|^2) where C is the number of edges in the graph
        Worst: O(|D| * |C|^2) where D is the number of data centres

    Space complexity:
        Input: O(|C| + |D|) where D is the number of data centres
        Aux: O(|C| + |D|)
    """

    verNum = len(maxIn)

    fN = FlowNetwork(verNum)
    fN.add_Connection(connections, maxIn, maxOut, targets)
    max_flow = fN.ford_fulkerson(origin, fN.imag_end_index)
    return max_flow




#Assignment2 Question2
class Node:
    """
    Written by LIANG DIZHEN
    Node class is to store the information of each node in the trie

    """
    def __init__(self, size = 27):
        """
            Written by LIANG DIZHEN
            Node class is to store the information of each node in the trie

            Precondition:
            Postcondition:

            Input:
                size - size of the link array (1 terminal node + 26 lowercase characters)
            Return:
                None

            Time complexity:
                Best: O(1)
                Worst: O(1)
            Space complexity:
                Input: O(1)
                Aux: O(1)
            """
        #terminal $ at index 0 +  lp + up
        self.link = [None] * size #link stores $ + 27 characters in order
        #frequncy of the terminal node
        self.freq = 0

        #to point to the terminal node that has the highest frequencies
        self.freqNode = None

        self.partialStr = None

        #previous node
        self.previous = None


class CatsTrie:
    def __init__(self, sentences):
        """
        Written by LIANG DIZHEN
        This __init__ function is to build a trie for the sentences list by inserting sentencce one by one


        Precondition:
        Postcondition:

        Input:
            sentences - sentences list
        Return:
            None
        Time complexity:
            Best: O(NM)
            Worst: O(NM) where N is the number of sentences and M is the length of the longest sentence
        Space complexity:
            Input: O(NM)
            Aux: O(NM)
        """
        self.root = Node()
        self.root.partialStr= ""

        #build trie
        for s in sentences:
            self.insert_recur(s)


    #recursion for each sentence in sentences
    def insert_recur(self, sentence):
        """
        Written by LIANG DIZHEN
        This insert_recur function is a starter of the actual recursion function that insert each char in the key
        into the trie.

        Precondition:
        Postcondition:

        Input:
            sentence - a sentence in the sentences list
        Return:
            None
        Time complexity:
            Best: O(M)
            Worst: O(M) where M is the length of the longest sentence
        Space complexity:
            Input: O(M)
            Aux: O(M)
        """
        current = self.root
        #lengh of sentence
        count = len(sentence)
        self.insert_recur_aux(current, sentence, count)
        #compare and record most frequent node from start
        #most frequent node vary during insertion given differnt prefix

        #for checking bug
        #print("\ninsert for " + key + " done\n")


    #recursion for each char in key
    def insert_recur_aux(self, current, key, count):
        """
        Written by LIANG DIZHEN
        Modified from the insert function in the lecture note, which is a recursion function that insert each char in the key
        into the trie and update the frequency of the terminal node and the most frequent node for each prefix node at the
        same time.

        Base Case: when count == 0, which means the end of the key is reached, if the sentence not yet has the terminal node,
        then create one terminal node, and if terminal node already exist both then update the frequency of the terminal node and the most frequent node for
        each prefix node at the same time. If the current terminal node has higher frequency than the most frequent node
        for the prefix node given the same prefix, then update the most frequent node for the prefix node to be the current
        terminal node. If the current terminal node has the same frequency as the most frequent node for the prefix node, then
        compare the string of the current terminal node and the most frequent node for the prefix node to check which one is
        lexicographically smaller, then update the most frequent node for the prefix node to be the current terminal node.
        If the current terminal node has lower frequency than the most frequent node for the prefix node, then do nothing.
        After that, update the most frequent node of the prefix node to be the current terminal node.

        Recursion Case: when count != 0, which means the end of the key is not reached, then check if the current char in the
        key is in the trie. If it is in the trie and there is next node, then go to the next node and continue the recursion with
        smaller key input. If the character not in the trie given the prefix then create a new node and concatenate the
        partial string from the previous node with the current character in key and assign to the new created node
        partial string from the previous node with the current character in key and assign to the new created node
        which is the next node of the current node which is pointed by the node'link list or array
        which has index (character ASCII value - 97 + 1) then go to the next node and continue the
        recursion with smaller key input.

        Precondition:
        Postcondition:

        Input:
            current - the current node
            key - the sentence to be inserted
            count - the character index in the key
        Return:
            None

        Time complexity:
            Best: O(NM)
            Worst: O(NM) where N is the number of sentences and M is the length of the longest sentence
        Space complexity:
            Input: O(NM)
            Aux: O(NM)
        """
        previous = current

        #count for the first created terminal node
        #base case O(1)
        if count == 0:
            # reach the end of the key
            #if the terminal node exist
            #update frequncy
            if current.link[0] is not None:
                current = current.link[0]
                #update frequency for the existing sentence at terminal node
                current.freq += 1
                #print(current.partialStr +  ": freq: " + str(current.freq))

                # update frequncy of current most frequent node for relavent prefix nodes
                #not None for it to go until root Node, but not update root node (inclusively)
                while current.previous is not None:
                    if current.freqNode.freq > current.previous.freqNode.freq:
                        current.previous.freqNode = current.freqNode
                        #print("current node: :" + current.partialStr + " freqNode: " + current.freqNode.partialStr)
                    elif current.freqNode.freq == current.previous.freqNode.freq:
                        #if the frequency the same, return the lexicographically smaller one
                        if current.freqNode.partialStr < current.previous.freqNode.partialStr:
                            current.previous.freqNode = current.freqNode
                    current = current.previous

            #when end of the key and no terminal node
            #create terminal node for every complete sentence
            
            else:
                index = 0
                current.link[index] = Node()
                current = current.link[index]
                current.partialStr = previous.partialStr

                #itself is always the most frequent node to itself
                #since given different prefix, itself is the most frequent node to itself
                current.freqNode = current
                #update previous node
                current.previous = previous
                #update frequncy
                current.freq += 1


                # if the previous node has no most frequent node since the sentence terminal node is firsly
                # created, for all the prefix node, compare the frequency of the current terminal node and the
                # most frequent node for the prefix node, if the frequency of the current terminal node is higher
                # update the freqNode of the prefix node to be the current terminal node, if the frequency of the
                # current terminal node is the same as the most frequent node for the prefix node, then compare the
                # string of the current terminal node and the most frequent node for the prefix node to check which
                # one is lexicographically smaller, then update freqNode to be the lecicographically smaller one.
                    #to run until root node
                while current.previous is not None:

                    #if the previous node has no freqNode
                    if current.previous.freqNode is None:
                        current.previous.freqNode = current.freqNode

                    else:
                        if current.freqNode.freq > current.previous.freqNode.freq:
                            current.previous.freqNode = current.freqNode
                            #print("current node: :" + current.partialStr + " freqNode: " + current.freqNode.partialStr)
                        elif current.freqNode.freq == current.previous.freqNode.freq:
                            #if the frequency the same, return the lexicographically smaller one
                            if current.freqNode.partialStr < current.previous.freqNode.partialStr:
                                current.previous.freqNode = current.freqNode

                    current = current.previous

                #print("terminalCount: " + str(terminalCount))



        #Recursive case starts from here !!!!!!!!!!!!!!!!!!!!!!!
        #O(n)
        #directly go to here as lone as there is char in key
        else:
            # find the index of the char with the ASCII value
            #since the link store them in order
            #alway get the first one for recursion
            #since every time input is 1 less
            #adjument for unchanged input with count
            #print(count)
            #key input keep unchanged, use count to adjust to correct first index of character
            adj_count = len(key) - count
            index = ord(key[0 + adj_count]) - 97 + 1

            #if the char node exist go alone the path
            if current.link[index] is not None:
                current = current.link[index]


                #the most frequent node change during the insertion
                #update the freqNode


            #if char node does not exist
            #create new node
            else:
                #print("New node created: " + current.partialStr)
                current.link[index] = Node()
                current = current.link[index]

                #store previous node
                current.previous = previous

                # store partial sentence
                current.partialStr = previous.partialStr + key[0 + adj_count]

            # recur with smaller input
            # using count to adjust the input
            count -= 1
            self.insert_recur_aux(current, key, count)




    #auto-complete
    def autoComplete(self, key):
        """
        Written by LIANG DIZHEN
        This function is to auto-complete the sentence with the given prefix
        If the sentence has not terminal node in the trie, return None
        If the sentence has terminal node in the trie, and the frequency of the sentence > frequency of the most frequent
        sentence with the same prefix, return the sentence
        If the sentence has terminal node in the trie, and the frequency of the sentence = frequency of the most frequent
        sentence with the same prefix, return the lexicographically smaller sentence
        If the sentence has terminal node in the trie, and the frequency of the sentence < frequency of the most frequent
        sentence with the same prefix, return the most frequent sentence with the same prefix
        If the sentence has no prefix in the trie, return most frequent sentence in the trie

        Precondition:
        Postcondition:

        Input:
            key: the sentence
        Return:
            a sentence with the given prompt by aviding the above conditions

        Time complexity:
            Best: O(X) given the sentence does not exist in the trie given the prompt, X is the length of the prompt
            Worst: O(X+Y) X is the length of the prompt, Y is the length of the most frequent sentence begin with the prompt
            as need to compare the strings' overall lexigraphical order
        Space complexity:
            Input: O(X) X is the length of the prompt
            Aux:  O(1)
        """
        # start ffrom the root
        current = self.root

        #go through the key 1 by 1
        for char in key:

            # calculate index
            # $ = 0, a = 1, b = 2, c = 3
            index = ord(char) - 97 + 1 #shift by 1

            #if path exit
            #
            if current.link[index] is not None:
                current = current.link[index]


            #if path does not exist
            #at lone as it has no way out, return None
            else:
                return None



        #move until to check when this sentence has terminal node or not
        #now we are at the leaf node(terminal)
        #print(current.level)

        # after finish the given string (list of chars)
        # check if there is a $
        # go through the terminal $, index = 0
        #read the end of the key
        # terminal $ at index 0
        tIndex = 0
        if current.link[tIndex] is not None:
            #if the freqNode has highest frequncy return it
            if current.freqNode.freq > current.link[tIndex].freq:
                return current.freqNode.partialStr
            
            # #if the frequency the same, return the lexicographically smaller one
            # if the prefix is the same, but not terminal node
            # return the most frequent node given the prefix
            elif current.freqNode.freq == current.link[tIndex].freq:
                #if the frequency the same, return the lexicographically smaller one
                if current.freqNode.partialStr < current.link[tIndex].partialStr:
                    return current.freqNode.partialStr
                else:
                    current = current.link[tIndex]
                    return current.partialStr
            
            else:
                #traverse to terminal node (leaf node) (index 0 of link)
                current = current.link[tIndex]
                return current.partialStr

            #return complete sentence from the terminal node
            #return current.partialStr
        
        #if it has not terminal node
        else:
            return current.freqNode.partialStr



# if __name__ == '__main__':
#     sentences = ['qgqlmiux', 'g', 'ygxikxhde', 'pszueamw', 'xfssuxcf', 'ng', 'uynv', 'juwjv', 'jt', 'qbpqflfjby',
#                  'qveqjbzr', 'igevwuat', 'qveqjbzr', 'mzsktobx', 'gouhmktouu', 'nlstgvxzqi', 'xfssuxcf', 'jfzaipjrzx',
#                  'sijcnfhhah', 'ohdhswcxuv', 'owxgnxtmp', 'g', 'yegoxnftr', 'dcktqp', 'ye', 'uqz', 'mzsktobx', 'mqiln',
#                  'xlxzqmtzh', 'pybplk', 'ng', 'evuw', 'qveqjbzr', 'jt', 'qveqjbzr', 'qstxov', 'swqkbd', 'ye',
#                  'ntgruagd', 'ohdhswcxuv', 'uf', 'mxxgmoyrwv', 'sijcnfhhah', 'fcbf', 'uynv', 'fcozq', 'uynv',
#                  'mxxgmoyrwv', 'uzdgivw', 'uqz', 'swqkbd', 'mxxgmoyrwv', 'wkltp', 'fcbf', 'bgayig', 'qgqlmiux',
#                  'dexasozo', 'pszueamw', 'fcozq', 'trmazkn', 'jfzaipjrzx', 'gouhmktouu', 'ng', 'depv', 'fcbf',
#                  'ygxikxhde', 'dexasozo', 'swqkbd', 'xfssuxcf', 'mqiln', 'ygxikxhde', 'n', 'wdvu', 'uynv', 'juwjv',
#                  'ccf', 'qstxov', 'wdvu', 'q', 'pybplk', 'wiiz', 'uh', 'yegoxnftr', 'lioky', 'n', 'nlstgvxzqi',
#                  'ibohwfyicj', 'yegoxnftr', 'hayzmd', 'trmazkn', 'n', 'dcktqp', 'uzdgivw', 'igevwuat', 'uh', 'uf',
#                  'mqiln', 'jt', 'swqkbd', 'xlxzqmtzh']
#     trie = CatsTrie(sentences)
#     t = trie.autoComplete("dck")

#     print(t)






















