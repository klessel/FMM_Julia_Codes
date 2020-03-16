function FMM_Construction2(tree::Tree,r,f::Function)
#This is a memory efficient, two pass algorithm that takes in a Tree
#which is a partition tree for a given matrix which is described by the
#function, f, below and returns its corresponding FMM representation.
#Computation of U,R,V,W is done via deepest first post-ordering.  B
#matrices are computed by in descending order (starting at root and finishing
#at the child).  Relavant matrices are then multiplied and added
#from child to root in order to compute each B matrix.  (Bottom Up Routine
# to compute B).
#INPUT:             r       (Int64) largest allowable rank of hankel blocks -
#                           this determines the amount of compression, and
#                           should be compatible with your input tree.
#                           (if the maximum allowable rank is p, then the
#                           tree partitions should be no smaller than 3p)
#                           Corresponding singular values below this rank
#                           will be dropped.
#                   tree    (Tree) contains partition dimensions for
#                           each subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#OUTPUT:            fmm     (FMM) contains matrices (Us, Vs, Bs Rs and
#                           Ws) that compose the fmm representation of the
#                           input matrix, in addition to the partition
#                           dimensions
#                   peakMem (Array{Int64,1}) first entry is a memory counter
#                           and second entry is the maximum memory
# Author: Kristen Lessel - November 2015
    N = tree.m;

    empty_tree = create_empty_tree(tree.depth);
    #pm_count = 0;
    #pm_max = 0;
    #peakMem = [pm_count,pm_max];
    fmm, dummy1, dummy2 = FMM_Basis_TranslationOp(empty_tree,tree,empty_tree,N,r,f);

    row_idx = 1;
    col_idx = 1;
    fmm  = FMM_Expansion_Coeffs_Diag(fmm,row_idx,col_idx,N,f);

    fmm
end


function FMM_Basis_TranslationOp(treeL::Tree,tree::Tree,treeR::Tree,N::Int64,r::Int64,f::Function)
#Computes U's, V's, R's, W's and D's of FMM structure, and stores these
#heirarchically
#INPUT:             N       (Int64) grid size
#                   rank    (Int64) largest allowable rank of hankel blocks -
#                           this determines the amount of compression, and
#                           should be compatible with your input tree.
#                           (if the maximum allowable rank is p, then the
#                           tree partitions should be no smaller than 3p)
#                           Corresponding singular values below this rank
#                           will be dropped.
#                   tree    (Union(Node_Spine,Leaf_Spine))
#                           contains partition dimensions for each
#                           subdivision of the input matrix, for each
#                           of which there are a left and right child
#                           structure
#OUTPUT:            fmm     (Union(Node,Leaf)) contains matrices
#                           (U, V, R and W) that compose the
#                           hss representation of the input matrix, in
#                           addition to the partition dimensions
#                   rH      (Array{Float64}(2)) row hankel block with 'current' U
#                           removed
#                   cH      (Array{Float64}(2)) column hankel block with 'current' V'
#                           removed
    if isa(treeL,Leaf_Spine) && isa(tree,Leaf_Spine) && isa(treeR,Leaf_Spine)

        d1_idx = collect(treeL.col_idx:treeL.col_idx+treeL.m-1);

        d2_idx = collect(tree.col_idx:tree.col_idx+tree.m-1);


        d3_idx = collect(treeR.col_idx:treeR.col_idx+treeR.m-1);

        m_idx = collect(treeL.col_idx:(treeL.col_idx+treeL.m+tree.m+treeR.m)-1);
        #generate indices for current off diagonal hankel block
        butm_idx = [1:m_idx[1]-1; m_idx[end]+1:N];

        diag_idx = collect(tree.col_idx:(tree.col_idx+tree.m)-1);

        #function call to generate lowest level row and column hankel blocks
        rH2 = f(diag_idx,butm_idx,N);
        cH2= f(butm_idx,diag_idx,N);

        #take svds of upper row/left column and lower row/right column hankel blocks
        Ur, Rr, Vr = svd_ranktol(rH2,r);
        Uc, Rc, Vc = svd_ranktol(cH2,r);

        rH2 = Array{Float64}(undef,0);
        cH2 = Array{Float64}(undef,0);

        leaf = Leaf(tree.m,tree.col_idx,tree.depth,Ur,Vc,f(d2_idx,d1_idx,N),
                    f(d2_idx,d2_idx,N),f(d2_idx,d3_idx,N));

        rH = diagm(Rr)*Vr';
        cH = Uc*diagm(Rc);

        leaf, rH, cH

    elseif isa(tree, Node_Spine)
        #If the center tree is a node

        #Descend into the tree, depth first
        if tree.treeL.depth >= tree.treeR.depth #left node deeper than right

            if isa(treeL,Node_Spine) && isa(tree,Node_Spine) && isa(treeR,Node_Spine)
                # 2) All nodes

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL.treeR,tree.treeL,tree.treeR,N,r,f);

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR.treeL,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeL.treeL.m + treeR.treeL.m + treeR.treeR.m;
                diag_width_r = treeL.treeL.m + treeL.treeR.m + treeR.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.treeR.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.treeL.m))];

            elseif isa(treeL,Node_Spine) &&isa(tree,Node_Spine) && isa(treeR,Leaf_Spine)
                # 3) Left and center are nodes, right is leaf
                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL.treeR,tree.treeL,tree.treeR,N,r,f);

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeL.treeL.m + treeR.m;
                diag_width_r = treeL.treeL.m + treeL.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.treeR.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.m))];

            elseif isa(treeL,Leaf_Spine) && isa(tree,Node_Spine) &&  isa(treeR,Leaf_Spine)
                # 8) Left and Right are leaves, center is node

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL,tree.treeL,tree.treeR,N,r,f);

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeR.m
                diag_width_r = treeL.m

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.m))];

            elseif isa(treeL,Leaf_Spine) && isa(tree,Node_Spine) && isa(treeR,Node_Spine)
                # 5) Left is a leaf, center and right are nodes

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL,tree.treeL,tree.treeR,N,r,f);

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR.treeL,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeR.treeL.m + treeR.treeR.m;
                diag_width_r = treeL.m + treeR.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.treeL.m))];

            end

        else #right node deeper than left
            if isa(treeL,Node_Spine) && isa(tree,Node_Spine) && isa(treeR,Node_Spine)
                # 2) All nodes

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR.treeL,N,r,f);

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL.treeR,tree.treeL,tree.treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeL.treeL.m + treeR.treeL.m + treeR.treeR.m;
                diag_width_r = treeL.treeL.m + treeL.treeR.m + treeR.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.treeR.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.treeL.m))];

            elseif isa(treeL,Node_Spine) && isa(tree,Node_Spine) && isa(treeR,Leaf_Spine)
                # 3) Left and center are nodes, right is leaf

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR,N,r,f);

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL.treeR,tree.treeL,tree.treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeL.treeL.m + treeR.m;
                diag_width_r = treeL.treeL.m + treeL.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.treeR.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.m))];

            elseif isa(treeL,Leaf_Spine) && isa(tree,Node_Spine) && isa(treeR,Leaf_Spine)
                # 8) Left and Right are leaves, center is node

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR,N,r,f);

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL,tree.treeL,tree.treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeR.m
                diag_width_r = treeL.m

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.m))];

            elseif isa(treeL,Leaf_Spine) && isa(tree, Node) && isa(treeR,Node_Spine)
                # 5) Left is a leaf, center and right are nodes

                #right recursive call
                fmmR, lowerRow, rightCol = FMM_Basis_TranslationOp(tree.treeL,tree.treeR,treeR.treeL,N,r,f);

                #left recursive call
                fmmL, upperRow, leftCol = FMM_Basis_TranslationOp(treeL,tree.treeL,tree.treeR,N,r,f);

                #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
                diag_width_l = treeR.treeL.m + treeR.treeR.m;
                diag_width_r = treeL.m + treeR.treeR.m;

                #indexes used to obtain compressed hankel blocks for the next level
                uR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_l):(N-(treeL.m + tree.treeL.m + tree.treeR.m))];
                lR_idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width_r):(N-(tree.treeL.m + tree.treeR.m + treeR.treeL.m))];

            end
        end

        rH_top = upperRow[:, uR_idx];
        rH_bottom = lowerRow[:, lR_idx];
        cH_left = leftCol[uR_idx,:];
        cH_right = rightCol[lR_idx,:];

        #merge remainder of Hankel blocks from children
        rH2 = [rH_top; rH_bottom];
        cH2 = [cH_left cH_right];

        #take svd of remaining portions of blocks to determine Rs and Ws
        Ur, Sr, Vr = svd_ranktol(rH2,r);
        Uc, Sc, Vc = svd_ranktol(cH2,r);

        rH2 = Array{Float64}(undef,0);
        cH2 = Array{Float64}(undef,0);

        #partition Us to get R's, partition V's to get Ws
        Rl = Ur[1:size(upperRow,1),:];
        Rr = Ur[size(upperRow,1)+1:end,:];
        Wl = Vc[1:size(leftCol,2),:];
        Wr = Vc[size(leftCol,2)+1:end,:];

        emptyleaf = Leaf(-1,-1,-1,Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0));
        fmm = Node(fmmL,fmmR,emptyleaf,emptyleaf,tree.m,tree.m,tree.col_idx,tree.col_idx,tree.depth,Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Array{Float64,2}(undef,0,0),Rl,Rr,Wl,Wr);

        #return relavant portions of hankel blocks
        rH = diagm(Sr)*Vr';
        cH = Uc*diagm(Sc);

        fmm, rH, cH

    elseif isa(tree, Leaf_Spine)
        #if center tree is a leaf, but left and/or right tree are/is node(s)

        if isa(treeL,Node_Spine) &&isa(tree,Leaf_Spine) && isa(treeR,Node_Spine)
            # 4) Left and Right trees are nodes, center tree is a leaf

            #recursive call
            fmm, row, col = FMM_Basis_TranslationOp(treeL.treeR,tree,treeR.treeL,N,r,f);

            # #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
            diag_width = treeL.treeL.m + treeR.treeR.m;

            # #indexes used to obtain compressed hankel blocks for the next level
            idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width):(N-(treeL.treeR.m + tree.m + treeR.treeL.m))];

        elseif isa(treeL,Leaf_Spine) &&isa(tree,Leaf_Spine) && isa(treeR,Node_Spine)
            #6) if left and center trees are leaves and right tree is a node

            #recursive call
            fmm, row, col = FMM_Basis_TranslationOp(treeL,tree,treeR.treeL,N,r,f);

            #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
            diag_width = treeR.treeR.m;

            #indexes used to obtain compressed hankel blocks for the next level
            idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width):(N-(treeL.m + tree.m + treeR.treeL.m))];

        elseif isa(treeL,Node_Spine) &&isa(tree,Leaf_Spine) && isa(treeR,Leaf_Spine)
            #7) if left tree is a node and center and right trees are leaves

            #recursive call
            fmm, row, col = FMM_Basis_TranslationOp(treeL.treeR,tree,treeR,N,r,f);

            #these are used to indexes of rol/col hankel blocks excluding corresponding part of diagonal block
            diag_width = treeL.treeL.m;

            #indexes used to obtain compressed hankel blocks for the next level
            idx = [1:(treeL.col_idx-1); (treeL.col_idx + diag_width):(N-(treeL.treeR.m +tree.m + treeR.m))];

        end

        rH = row[:,idx];
        cH = col[idx,:];

        fmm, rH, cH
end
end

#function FMM_Expansion_Coeffs_Diag(fmm::Fss,row_idx::Int64,col_idx::Int64,N::Int64,peakMem::Array{Int64,1})
function FMM_Expansion_Coeffs_Diag(fmm::FMM,row_idx::Int64,col_idx::Int64,N::Int64,f::Function)
#Recursively computes Expansion Coefficients (B's) of FMM structure for the
#upper left and lower right block of the current node
#
#   ---------------
#  \   \   \B13\B14\
#  \-------\-------\
#  \   \   \   \B24\
#  \---------------\
#  \B31\   \   \   \
#  \-------\-------\
#  \B41\B42\   \   \
#   ---------------
#INPUT:             fmm     (Union(Node,Leaf)) branch of current node,
#                   row_idx (Int64) row index of current node
#                   col_idx (Int64) col index of current node
#                   N       (Int64) grid size
#OUTPUT:            fmm     (Union(Node,Leaf)) contains matrices (U, V, B R and
#                           W) that compose the HSS representation of the
#                           input tree for a corresponding matrix, as well
#                           as corresponding partition dimensions and
#                           depth for each node

    if isa(fmm.fmmUL,Leaf) && isa(fmm.fmmLR,Leaf)
        #do nothing.  Return
        #   -------
        #  \ D \ D \
        #  \-------\
        #  \ D \ D \
        #   -------

        fmm

    elseif !isa(fmm.fmmUL,Leaf) && !isa(fmm.fmmLR,Leaf)

        #Off diagonal block calls.  Recursively compute B13,B14,B24,B31,B41,B42
        #   -------
        #  \   \ x \
        #  \-------\
        #  \ x \   \
        #   -------
        col_idx = col_idx + fmm.fmmUL.m;
        fmmUR, fmmLL = FMM_Expansion_Coeffs_OffDiag(fmm.fmmUL,fmm.fmmLR,row_idx,col_idx,N,f);
        col_idx = col_idx - fmm.fmmUL.m;

        #Upper Left Diagonal Block FMM Call
        #   -------
        #  \ x \   \
        #  \-------\
        #  \   \   \
        #   -------
        fmmUL = FMM_Expansion_Coeffs_Diag(fmm.fmmUL,row_idx,col_idx,N,f);

        #Lower Right Diagonal Block FMM Call
        #   -------
        #  \   \   \
        #  \-------\
        #  \   \ x \
        #   -------
        col_idx = col_idx + fmm.fmmUL.m;
        row_idx = row_idx + fmm.fmmUL.m;
        fmmLR = FMM_Expansion_Coeffs_Diag(fmm.fmmLR,row_idx,col_idx,N,f) ;

        fmm.fmmUR = fmmUR;
        fmm.fmmLL = fmmLL;
        fmm.fmmUL = fmmUL;
        fmm.fmmLR = fmmLR;
        fmm

    elseif isa(fmm.fmmUL,Leaf) && !isa(fmm.fmmLR,Leaf)

        #compute both upper b's and lower b's here (2 of them)
        col_idx = col_idx + fmm.fmmUL.m;
        fmmUR, fmmLL = FMM_Expansion_Coeffs_OffDiag(fmm.fmmUL,fmm.fmmLR,row_idx,col_idx,N,f);
        col_idx = col_idx - fmm.fmmUL.m;

        #Upper Left Diagonal BLock is a Leaf, so no call
        #fmmUL  = LeafOffDiag(fmm.fmmUL.m,fmm.fmmLR.m,fmm.fmmUL.col_idx,fmm.fmmLR.col_idx,0,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));

        #Lower Right Diagonal Block FMM Call
        col_idx = col_idx + fmm.fmmUL.m;
        row_idx = row_idx + fmm.fmmUL.m;
        fmmLR = FMM_Expansion_Coeffs_Diag(fmm.fmmLR,row_idx,col_idx,N,f);

        fmm.fmmUR = fmmUR;
        fmm.fmmLL = fmmLL;

        fmm.fmmLR = fmmLR
        fmm

    elseif !isa(fmm.fmmUL,Leaf) && isa(fmm.fmmLR,Leaf)

        #compute both upper b's and lower b's here (2 of them)
        col_idx = col_idx + fmm.fmmUL.m;
        fmmUR, fmmLL = FMM_Expansion_Coeffs_OffDiag(fmm.fmmUL,fmm.fmmLR,row_idx,col_idx,N,f);
        col_idx = col_idx - fmm.fmmUL.m;

        #Upper Left Diagonal Block FMM Call
        fmmUL = FMM_Expansion_Coeffs_Diag(fmm.fmmUL,row_idx,col_idx,N,f);

        #Lower Right Diagonal Block is a Leaf, so no call
        #fmmLR  = LeafOffDiag(fmm.fmmLR.m,fmm.fmmUL.m,fmm.fmmLR.col_idx,fmm.fmmUL.col_idx,0,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));


        fmm.fmmUR = fmmUR;
        fmm.fmmLL = fmmLL;
        fmm.fmmUL = fmmUL;

        fmm

    end

end

function FMM_Expansion_Coeffs_OffDiag(fmmUL::FMM,fmmLR::FMM,row_idx::Int64,col_idx::Int64,N::Int64,f::Function)

    if isa(fmmUL,Leaf) && isa(fmmLR,Leaf)
        ## These are zero arrays becuase B's at the leaf will be computed by Off_Diag_Outter.

        leafUR = LeafOffDiag(fmmUL.m,fmmLR.m,fmmUL.col_idx,fmmLR.col_idx,0,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));

        leafLR = LeafOffDiag(fmmLR.m,fmmUL.m,fmmLR.col_idx,fmmUL.col_idx,0,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));

        leafUR, leafLR

    elseif !isa(fmmUL,Leaf) && !isa(fmmLR,Leaf)
        #upper left corner of both upper right and lower left blocks
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \       \ x \<--\-- MxN
        #   \        -------\
        #   \       \   \   \
        #   \---------------\
        #   \ x \<--\-------\-- NxM
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (1)


        B13, B31 = FMM_Expansion_Coeffs_OffDiag_Outter(fmmUL.fmmUL,fmmLR.fmmUL,row_idx,col_idx,N,f);

        #Lower left corner lower left block, and also the upper right corner of
        #the upper right block
        #    ---------------
        #   \       \   \ x \
        #   \       \---\---\
        #   \       \   \   \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \ x \   \       \
        #    ---------------
        #           (2)

        col_idx = col_idx + fmmLR.fmmUL.m;
        B14, B41 = FMM_Expansion_Coeffs_OffDiag_Outter(fmmUL.fmmUL,fmmLR.fmmLR,row_idx,col_idx,N,f);

        #lower right corner of both lower left and upper left blocks
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \   \ x \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \   \ x \       \
        #    ---------------
        #           (3)

        row_idx = row_idx +fmmUL.fmmUL.m;
        B24, B42 = FMM_Expansion_Coeffs_OffDiag_Outter(fmmUL.fmmLR,fmmLR.fmmLR,row_idx,col_idx,N,f);

        #upper right corner of lower block, lower left corner of upper block
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \ x \   \
        #   \---------------\
        #   \   \ o \       \
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (4)

        col_idx = col_idx - fmmLR.fmmUL.m #changed on 4/18/16
        fmmUR_LL, fmmLL_UR = FMM_Expansion_Coeffs_OffDiag(fmmUL.fmmLR,fmmLR.fmmUL,row_idx,col_idx,N,f); # x, o =  FMM_Expansion_Coeffs_OffDiag()

        offDiag_depth = max(fmmLR.depth,fmmUL.depth);
        emptyleaf = LeafOffDiag(-1,-1,-1,-1,-1,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));
        fmmUR = Node(emptyleaf, emptyleaf, emptyleaf,fmmUR_LL,fmmUL.m,fmmLR.m,fmmLR.col_idx,
                     fmmUL.col_idx,offDiag_depth ,B13,B14,B24,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0));
        fmmLL = Node(emptyleaf,emptyleaf,fmmLL_UR,emptyleaf,fmmLR.m,fmmUL.m,fmmUL.col_idx,
                     fmmLR.col_idx,offDiag_depth, Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     B31,B41,B42,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0));

        fmmUR, fmmLL

    elseif isa(fmmUL,Leaf) && !isa(fmmLR,Leaf)
        #right block of upper right corner and bottom block of lower left corner
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \       \   \   \
        #   \   D   \   \B13\<-- MxN
        #   \       \   \   \
        #   \---------------\
        #   \       \   \   \
        #   \-------\---\---\
        #   \  B31  \   \   \
        #    ---------------
        #       ^    (5)
        #       \
        #      NxM

        col_idx = col_idx +fmmLR.fmmUL.m;
        B13, B31 = FMM_Expansion_Coeffs_OffDiag_Outter(fmmUL,fmmLR.fmmLR,row_idx,col_idx,N,f);

        #left block of upper right corner and top block of lower left corner
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \       \   \   \
        #   \   D   \ x \   \<-- MxN
        #   \       \   \   \
        #   \---------------\
        #   \   x   \   \   \
        #   \-------\---\---\
        #   \       \   \   \
        #    ---------------
        #       ^    (6)
        #       \
        #      NxM


        col_idx = col_idx - fmmLR.fmmUL.m;

        fmmUR_L, fmmLL_L = FMM_Expansion_Coeffs_OffDiag(fmmUL,fmmLR.fmmUL,row_idx,col_idx,N,f);
        emptyleaf = LeafOffDiag(-1,-1,-1,-1,-1,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));
        fmmUR = Node(emptyleaf, emptyleaf, emptyleaf,fmmUR_L,fmmUL.m,fmmLR.m,fmmLR.col_idx,
                     fmmUL.col_idx,fmmLR.depth ,B13,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));
        fmmLL = Node(emptyleaf,emptyleaf,fmmLL_L,emptyleaf,fmmLR.m,fmmUL.m,fmmUL.col_idx,
                     fmmLR.col_idx,fmmLR.depth, Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     B31,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));

        fmmUR, fmmLL

    elseif !isa(fmmUL,Leaf) && isa(fmmLR,Leaf)
        #upper block of upper right corner, and left block of lower left corner
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \   \   \  B13  \
        #   \---\---\-------\
        #   \   \   \       \
        #   \---------------\
        #   \   \   \       \
        #   \B31\   \       \
        #   \   \   \       \
        #    ---------------
        #           (7)

        B13, B31 = FMM_Expansion_Coeffs_OffDiag_Outter(fmmUL.fmmUL,fmmLR,row_idx,col_idx,N,f);

        #left block of upper right corner and top block of lower left corner
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \   \   \       \
        #   \---\---\-------\
        #   \   \   \   x   \
        #   \---------------\
        #   \   \   \       \
        #   \   \ x \       \
        #   \   \   \       \
        #    ---------------
        #          (8)


        row_idx = row_idx + fmmUL.fmmUL.m;
        fmmUR_R, fmmLL_R = FMM_Expansion_Coeffs_OffDiag(fmmUL.fmmLR,fmmLR,row_idx,col_idx,N,f);

        emptyleaf = LeafOffDiag(-1,-1,-1,-1,-1,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));
        fmmUR = Node(emptyleaf, emptyleaf, emptyleaf,fmmUR_R,fmmUL.m,fmmLR.m,fmmLR.col_idx,
                     fmmUL.col_idx,fmmUL.depth ,B13,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));
        fmmLL = Node(emptyleaf,emptyleaf,fmmLL_R,emptyleaf,fmmLR.m,fmmUL.m,fmmUL.col_idx,
                     fmmLR.col_idx,fmmUL.depth, Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     B31,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                     Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0));

        fmmUR, fmmLL

    end
end

function FMM_Expansion_Coeffs_OffDiag_Outter(fmmL::FMM,fmmR::FMM,row_idx::Int64,col_idx::Int64,N::Int64,f::Function)
#function FMM_Expansion_Coeffs_OffDiag(treeL::Hss,treeR::Hss,row_idx::Int64,col_idx::Int64,N::Int64)
#Computes Expansion Coefficients (B's) of HSS structure for the upper right
#and lower left block of the current node
#   -------
#  \   \ x \
#  \-------\
#  \ x \   \
#   -------
#INPUT:             fmmL   (Union(Node,Leaf)) left branch of current node
#                   fmmR   (Union(Node,Leaf)) right branch of current node
#                   row_idx (Int64) row index of current node
#                   col_idx (Int64) col index of current node
#                   N       (Int64) grid size
#OUTPUT:            Bu      (Array{Float62,2}) Expansion coefficient matrix B at the
#                           current level which corresponds to the upper
#                           right block
#                   Bl      (Array{Float62,2}) Expansion coefficient matrix B at the
#                           current level which corresponds to the lower
#                           left block

    Au = Array{Float64,2};
    Bu = Array{Float64,2};
    Al = Array{Float64,2};
    Bl = Array{Float64,2};

    bu_l = Array{Float64,2};
    bu_r = Array{Float64,2};
    bl_l = Array{Float64,2};
    bl_r = Array{Float64,2};
    if isa(fmmL, Leaf) && isa(fmmR, Leaf) #If node is a leaf
        m_idx = collect((row_idx):(row_idx+fmmL.m-1));
        butm_idx = collect(col_idx:(col_idx+fmmR.m-1));

        Au = f(m_idx,butm_idx,N);
        Bu = fmmL.U'*Au*fmmR.V;

        Au = Array{Float64}(undef,0);

        Al = f(butm_idx,m_idx,N);
        Bl = fmmR.U'*Al*fmmL.V;

        Al = Array{Float64}(undef,0);

        Bu, Bl

    elseif !isa(fmmL, Leaf) && !isa(fmmR,Leaf)  #If neither node is a leaf
        #COMPUTE EXPANSION COEFFICIENTS

        #upper left corner of both upper right and lower left blocks
        #(these are computed at the same time. One block will be of
        #dimension MxN, the other will always be of dimension NxM)
        #
        #    ---------------
        #   \       \ x \<--\-- MxN
        #   \        -------\
        #   \       \   \   \
        #   \---------------\
        #   \ x \<--\-------\-- NxM
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (1)

        Bu_ll, Bl_ll = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmUL,fmmR.fmmUL,row_idx,col_idx,N,f);
        bu_ll = fmmL.Rl'*Bu_ll*fmmR.Wl;
        bl_ll = fmmR.Rl'*Bl_ll*fmmL.Wl;

        Bu_ll = Array{Float64}(undef,0);
        Bl_ll = Array{Float64}(undef,0);

        #Lower left corner lower left block, and also the upper right corner of
        #the upper right block
        #    ---------------
        #   \       \   \ x \
        #   \       \---\---\
        #   \       \   \   \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \ x \   \       \
        #    ---------------
        #           (2)

        col_idx = col_idx +fmmR.fmmUL.m;
        Bu_lr, Bl_lr = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmUL,fmmR.fmmLR,row_idx,col_idx,N,f);
        bu_lr = fmmL.Rl'*Bu_lr*fmmR.Wr;
        bl_lr = fmmR.Rr'*Bl_lr*fmmL.Wl;

        Bu_lr = Array{Float64}(undef,0);
        Bl_lr = Array{Float64}(undef,0);

        #upper right corner of lower block, lower left corner of upper block
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \ x \   \
        #   \---------------\
        #   \   \ x \       \
        #   \-------\       \
        #   \   \   \       \
        #    ---------------
        #           (3)

        col_idx = col_idx -fmmR.fmmUL.m;
        row_idx = row_idx +fmmL.fmmUL.m;
        Bu_rl, Bl_rl = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmLR,fmmR.fmmUL,row_idx,col_idx,N,f);
        bu_rl = fmmL.Rr'*Bu_rl*fmmR.Wl;
        bl_rl = fmmR.Rl'*Bl_rl*fmmL.Wr;

        Bu_rl = Array{Float64}(undef,0);
        Bl_rl = Array{Float64}(undef,0);

        #lower right corner of both lower left and upper left blocks
        #    ---------------
        #   \       \   \   \
        #   \       \---\---\
        #   \       \   \ x \
        #   \---------------\
        #   \   \   \       \
        #   \-------\       \
        #   \   \ x \       \
        #    ---------------
        #           (4)

        col_idx = col_idx +fmmR.fmmUL.m;
        Bu_rr, Bl_rr = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmLR,fmmR.fmmLR,row_idx,col_idx,N,f);
        bu_rr = fmmL.Rr'*Bu_rr*fmmR.Wr;
        bl_rr = fmmR.Rr'*Bl_rr*fmmL.Wr;

        Bu_rr = Array{Float64}(undef,0);
        Bl_rr = Array{Float64}(undef,0);

        #Add b's to form the B's
        Bu = bu_ll+bu_rl+bu_lr+bu_rr;
        Bl = bl_ll+bl_rl+bl_lr+bl_rr;

        bu_ll = Array{Float64}(undef,0);
        bu_rl = Array{Float64}(undef,0);
        bu_lr = Array{Float64}(undef,0);
        bu_rr = Array{Float64}(undef,0);
        bl_ll = Array{Float64}(undef,0);
        bl_rl = Array{Float64}(undef,0);
        bl_lr = Array{Float64}(undef,0);
        bl_rr = Array{Float64}(undef,0);

        Bu, Bl

    elseif isa(fmmL, Leaf) && !isa(fmmR, Leaf) #If left node is a leaf, and right node is not
        #COMPUTE EXPANSION COEFFICIENTS

        #left block of upper right corner and upper block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   \   \
        #   \       \ x \<--\-- MxN
        #   \       \   \   \
        #   \---------------\
        #   \   x<--\-------\-- NxM
        #   \-------\       \
        #   \       \       \
        #    ---------------
        #           (5)

        Bu_l, Bl_l = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL,fmmR.fmmUL,row_idx,col_idx,N,f);
        bu_l = Bu_l*fmmR.Wl;
        bl_l = fmmR.Rl'*Bl_l;

        Bu_l = Array{Float64}(undef,0);
        Bl_l = Array{Float64}(undef,0);

        #right block of upper right corner, upper left block of lower left corner.
        #    ---------------
        #   \       \   \   \
        #   \       \   \ x \
        #   \       \   \   \
        #   \---------------\
        #   \       \       \
        #   \-------\       \
        #   \   x   \       \
        #    ---------------
        #           (6)

        col_idx = col_idx + fmmR.fmmUL.m;
        Bu_r, Bl_r = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL,fmmR.fmmLR,row_idx,col_idx,N,f);
        col_idx = col_idx - fmmR.fmmUL.m;
        bu_r = Bu_r*fmmR.Wr;
        bl_r = fmmR.Rr'*Bl_r;

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        Bu = bu_l + bu_r;
        Bl = bl_l + bl_r;

        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);
        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);

        Bu, Bl

    elseif !isa(fmmL, Leaf) && isa(fmmR, Leaf) #If right node is a leaf, and left node is not
        #COMPUTE EXPANSION COEFFICIENTS

        #upper block of upper right corner and left block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \   x <-\-- MxN
        #   \       \-------\
        #   \       \       \
        #   \---------------\
        #   \   \   \       \
        #   \ x \<--\-------\-- NxM
        #   \   \   \       \
        #    ---------------
        #           (7)

        Bu_l, Bl_l = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmUL,fmmR,row_idx,col_idx,N,f);
        bu_l = fmmL.Rl'*Bu_l;
        bl_l = Bl_l*fmmL.Wl;

        Bu_l = Array{Float64}(undef,0);
        Bl_l = Array{Float64}(undef,0);

        #upper block of upper right corner and left block of lower left corner
        #(these are computed in the same pass)
        #    ---------------
        #   \       \       \
        #   \       \-------\
        #   \       \   x   \
        #   \---------------\
        #   \   \   \       \
        #   \   \ x \       \
        #   \   \   \       \
        #    ---------------
        #           (8)

        row_idx = row_idx + fmmL.fmmUL.m;
        Bu_r, Bl_r = FMM_Expansion_Coeffs_OffDiag_Outter(fmmL.fmmLR,fmmR,row_idx,col_idx,N,f);
        row_idx = row_idx - fmmL.fmmUL.m;
        bu_r = fmmL.Rr'*Bu_r;
        bl_r = Bl_r*fmmL.Wr;

        Bu_r = Array{Float64}(undef,0);
        Bl_r = Array{Float64}(undef,0);

        Bu = bu_l + bu_r;
        Bl = bl_l + bl_r;

        bu_l = Array{Float64}(undef,0);
        bu_r = Array{Float64}(undef,0);
        bl_l = Array{Float64}(undef,0);
        bl_r = Array{Float64}(undef,0);

        Bu, Bl
    end
end
