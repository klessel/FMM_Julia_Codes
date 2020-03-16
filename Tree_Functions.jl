function Gen_Tree_Complete(minPart::Int64,N::Int64,d::Int64)
#function(minPart,N)
#This code is for testing purposes only.  Generates a tree
#structure that is not associated with any function.  Function takes in the
#size of the desired matrix and generates a tree by splitting this
#repeatedly in half on the left and right for a 'complete' tree.
#INPUT:             minPart     (Int64) minimum number of partitions
#                               type of tree that will be generated
#                   N           (Uint) number of points in desired
#                               function
#OUTPUT:            tree        (Union(Node_Spine,Leaf_Spine) contains split dimensions and depth at
#                               each level

    if N/2 < minPart # leaf
        depth = 0;
        leafL = Leaf_Spine(N,d,depth);

        d = d + N;
        leafL, depth, d
    else # node
        treeL, depth_l, d_l = Gen_Tree_Complete(minPart,convert(Int64,N/2),d);
        treeR, depth_r, d_r = Gen_Tree_Complete(minPart,convert(Int64,N/2),d_l);
        
        depth = 1 + max(depth_l,depth_r);
        tree = Node_Spine(treeL,treeR,N,d,depth);
        tree, depth, d_r
    end
end

function Gen_Tree_Mtree(minPart::Int64,N::Int64)
#function [tree] = Gen_TestTree(tree,N)
#This code is for testing purposes only. Artificually generates a tree 
#structure that is not associated with any function.  Function takes in the
#number of desired nodes and generates a Worst case memory tree
# must have the line 'using Polynomials' in the main routine
#
#INPUT:             
#                   N           (Int64) number of points in desired 
#                               function
#                   minPart     (Int64) minimum number of partitions 
#OUTPUT:            tree        (Tree) contains split dimensions at   
#                               each level, as well as whether or not the 
#                               node is a leaf

    
    # number of nodes
    n = 2*floor(N/minPart) -1; # n = N_N = N_L -1, and N_L = N/minPart. 

    #depth of tree
    nodeRoots = roots(Poly([1-n,1,1]));
    d_exact  = nodeRoots[nodeRoots.>0]; #choose the positive root.
    
    #find maximum depth, d_max, of the tree we will generate
    d_max = convert(Int64,ceil(d_exact[1]));

    #Call the function initially on the root node of the tree.
    #   .<-
    #  / \
    # /\ /\ 
    #/\ /\/\
    tree = Gen_MTree_Right(N,minPart,d_max);
    tree
end

function Gen_MTree_Right(m,minPart,depth)
#Descend into right branch of the root of the current subtree (pictured below)   
#   .
#  / \<-
# /\ /\ 
#/\ /\/\
#INPUT:             
#                   m           (Int64) partition dimension of current  
#                               node
#                   minPart     (Int64) minimum number of partitions
#                   depth       (Int64) depth of current node
#                   pB          (Bool) a flag that denotes whether the
#                               current node is on the prime branch
#OUTPUT:            tree        (Tree) contains split dimensions at   
#                               each level, as well as whether or not the 
#                               node is a leaf
    if !(m < 2*minPart)
    #Case 1: If we do have enough rows to split into 2 blocks of minimum
    #partition size split and recurse on both children.  If we do not, then do
    #not partition further; Return two leaf nodes.

        #max number of blocks of minimum partition size we can have on this left subtree
        numBlocks = convert(Int64,floor(m/minPart)); 
        #remaining depth of the left subtree.
        r_depth = numBlocks-1; 

        #if we cannot generate a full left subtree - only generate children
        #corresponding to the the number of partitions we have left.
        if m <(depth +1)*minPart
            mL = r_depth*minPart;
        else
            mL = depth*minPart;
        end

        #Left Subtree Call
        treeL = Gen_MTree_Left(mL,minPart,depth);
        
        #Right Subtree Call
        #do we have enough rows to partition into two blocks of minumum partition size?
        if  m >= 2*minPart && m < 3*minPart
            #Case A: do we only have enough to split into two blocks and no 
            #more?(blocks must be of at least minimum partition size and no 
            #more than 2 times the minimum partition size)
            #Ex: minPart = 60
            #      /\170
            #   110  60
            #

            mL = m- minPart; #extra rows tacked onto left child
            leafL = Leaf_Spine(mL,-1,-1);

            mR = minPart;
            leafR = Leaf_Spine(mR,-1,-1);

            tree  = Node_Spine(leafL,leafR,m,-1,-1);
            tree
        elseif m<(depth+1)*minPart
            #Case B: do we have enough to split into more than two blocks, but cannot
            #generate a full left going subtree?
            #Ex:(N = 4096) Depth of full tree = 12. MinPart = 60. 
            #(Though here I am only showing partition  dimensions for 3 levels.
            #Dots indicate a part of the tree not shown. )
            #   . . /\
            #  .     /\204
            #  . 144/\ 60
            #  .  60  84
            #  /\

            mR = m - (numBlocks-1)*minPart;

            treeR = Leaf_Spine(mR,-1,-1);
            tree = Node_Spine(treeL,treeR,m,-1,-1);
            tree
        else #(depth+1)*minPart < m
            #Case C: We have enough rows to generate a full leftgoing subtree of depth d_max
            
            mR = m- depth*minPart;
            treeR = Gen_MTree_Right(mR,minPart,depth);

            tree = Node_Spine(treeL,treeR,m,-1,-1);
            tree
        end
        
    else #Case 2: we did not have enough rows to partition - label node as a leaf and return
        leaf = Leaf_Spine(m,-1,-1); 
        leaf
        
    end
    
end

function Gen_MTree_Left(m,minPart,depth)
#Generate left branch of the root of the current subtree (pictured below)    
#     .
#    / \
# ->/\ /\ 
#  /\ /\/\
#INPUT:             
#                   m           (Int64) partition dimension of current  
#                               node
#                   minPart     (Int64) minimum number of partitions
#                   depth       (Int64) depth of current node
#                   pB          (Bool) a flag that denotes whether the
#                               current node is on the prime branch
#OUTPUT:            tree        (Tree) contains split dimensions at   
#                               each level, as well as whether or not the 
#                               node is a leaf

    if m < 2*minPart  # if Leaf
       
        #pB = false;
        leafL = Leaf_Spine(m,-1,-1);

        leafL
    else #node
        depth -= 1;
        mL = m - minPart;
        treeL = Gen_MTree_Left(mL,minPart,depth);

        treeR = Leaf_Spine(minPart,-1,-1);

        tree = Node_Spine(treeL,treeR,m,-1,-1); 
        tree
    end
    
end

function label_tree(tree::Tree,col_idx)
#This function labels each node with the starting index (row/col value)
# of its corresponding diagonal block. Can also label the depth if the
# commented lines are uncommented
#
#INPUT:             tree    (Tree) contains partition dimensions for 
#                           each subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#                   col_idx       (Int64) row and col index which corresponds to 
#                           the first element in the current diagonal block 
#                           at every node.
#OUTPUT:            tree    (Tree) same as input,  that contains valid
#                           diagonal (row/col) index at each node
#                           (and depth if commented lines are uncommented.
#                           
# Author: Kristen Lessel - Sept 2014

    if !isa(tree,Leaf_Spine) # if not a leaf node
        
        tree.col_idx = col_idx;
        tree.treeL,col_idx = label_tree(tree.treeL,col_idx);
        depth_l = tree.treeL.depth;
        
        tree.treeR,col_idx = label_tree(tree.treeR,col_idx);
        depth_r = tree.treeR.depth;
        
        tree.depth = 1 + max(depth_l,depth_r);
        tree, col_idx
    else 
        tree.depth = 0;
        tree.col_idx = col_idx;
        col_idx = col_idx + tree.m;
        tree, col_idx
    end
end

function create_empty_tree(depth)
#creates a complete tree with partition dimensions of 0 at each node
#INPUT:                depth (Int64) denoting depth of tree)
#OUTPUT                empty_tree (Tree) empty partition tree with given depth
    if depth == 0
    empty_leaf = Leaf_Spine(0,1,depth)
    
    empty_leaf
    else
        depth = depth-1;
        empty_treeL = create_empty_tree(depth);
        empty_treeR = create_empty_tree(depth);

        empty_tree = Node_Spine(empty_treeL,empty_treeR,0,1,depth+1);
        empty_tree
    end
    
end

function create_empty_fmm(depth)
    #creates a complete fmm representation  with partition dimensions of 0 at each node
    #The fmmUR and fmmUL nodes are flipped for use with the FMM_Multiply()
    #Creates an empty FMM Representation of this form:
    #    ---------------
    #   \       \   \ x \
    #   \       \-------\  
    #   \       \   \   \
    #   \---------------\
    #   \   \   \       \
    #   \-------\       \
    #   \ x \   \       \
    #    ---------------
    #         
    #The nodes on the diagonal could effectively be empty because they aren't used,
    #but this code assigns empty nodes/leaves/arrays to thesse as well.
    #
    #    
    #INPUT:                depth (Int64) denoting depth of tree)
    #OUTPUT                empty_tree (Tree) empty partition tree with given depth
    
    if depth == 0
        empty_leaf = Leaf(-1,-1,depth,Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0));
        empty_leaf

    else
        depth = depth-1;
        empty_treeUR,empty_treeLL = create_empty_fmm_offdiag(depth);
        empty_treeL = create_empty_fmm(depth);
        empty_treeR = create_empty_fmm(depth);
        
        empty_tree = Node(empty_treeL,empty_treeR,empty_treeLL,empty_treeUR,-1,-1,-1,-1,depth+1,
                          Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                          Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                          Array(Float64,0,0),Array(Float64,0,0));
        
        empty_tree
    end
    
end

function create_empty_fmm_offdiag(depth)
    if depth == 0
        empty_leaf_offdiag = LeafOffDiag(-1,-1,-1,-1,0,Array(Float64,0,0),Array(Float64,0,0),
                                          Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0));
        
        empty_leaf_offdiag, empty_leaf_offdiag

    else

        depth = depth-1;
        empty_treeUR, empty_treeLL = create_empty_fmm_offdiag(depth);

        empty_leaf = LeafOffDiag(-1,-1,-1,-1,-1,Array(Float64,0,0),Array(Float64,0,0),
                                 Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0));
        empty_treeUR2 = Node(empty_leaf, empty_leaf, empty_leaf,empty_treeUR,-1,-1,-1,
                     -1,depth+1,Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                     Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                     Array(Float64,0,0),Array(Float64,0,0));
        empty_treeLL2 = Node(empty_leaf, empty_leaf,empty_treeLL,empty_leaf,-1,-1,-1,
                     -1,depth+1,Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                     Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),Array(Float64,0,0),
                     Array(Float64,0,0),Array(Float64,0,0));

        empty_treeUR2, empty_treeLL2
    end

end

function create_empty_gtree(depth)
#creates a complete tree with partition dimensions of 0 at each node
#INPUT:                depth (Int64) denoting depth of tree)
#OUTPUT                empty_tree (Tree) empty partition tree with given depth
    if depth == 0
        g = Array(Float64,0,1); #i have this as 3, but do i need to give an input fmm structure in general to match the size of the g's? Does the size of the empty matrix matter, is it used at all? 
        empty_gleaf = Leaf_Gtree(0,1,depth,g)
        
        empty_gleaf
    else
        depth = depth-1;
        empty_gtreeL = create_empty_gtree(depth);
        empty_gtreeR = create_empty_gtree(depth);

        g = Array(Float64,0,1);
        empty_gtree = Node_Gtree(empty_gtreeL,empty_gtreeR,0,1,depth+1,g);
        empty_gtree
    end
    
end

