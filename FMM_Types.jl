## FMM_Types.jl
# 1-16-18

module FMM_Types

export Leaf, LeafOffDiag, Node, FMM, Leaf_Spine, Node_Spine, Tree, Leaf_Gtree, Node_Gtree, Gtree

mutable struct Leaf{T1 <: Integer, T2 <: Number}
    m::T1 #dimension of dense diagonal (square) matrix at current node
    col_idx::T1 #column index number of the matrix at the current node
    depth::T1
    U::Array{T2,2}
    V::Array{T2,2}
    D_l::Array{T2,2}
    D_0::Array{T2,2}
    D_r::Array{T2,2}
end

mutable struct LeafOffDiag{T1 <: Integer, T2 <: Number}
    m::T1
    n::T1
    col_idx::T1
    row_idx::T1
    depth::T1
    B13::Array{T2,2}
    B14::Array{T2,2}
    B24::Array{T2,2}
    B31::Array{T2,2}
    B41::Array{T2,2}
    B42::Array{T2,2}
end

mutable struct Node{T1 <: Integer, T2 <: Number}
    fmmUL::Union{Node,Leaf,LeafOffDiag}
    fmmLR::Union{Node,Leaf,LeafOffDiag}
    fmmUR::Union{Node,Leaf,LeafOffDiag}
    fmmLL::Union{Node,Leaf,LeafOffDiag}
    m::T1
    n::T1
    col_idx::T1
    row_idx::T1
    depth::T1
    B13::Array{T2,2}
    B14::Array{T2,2}
    B24::Array{T2,2}
    B31::Array{T2,2}
    B41::Array{T2,2}
    B42::Array{T2,2}
    Rl::Array{T2,2}
    Rr::Array{T2,2}
    Wl::Array{T2,2}
    Wr::Array{T2,2}
end

#Define FMM structure
FMM = Union{Leaf,Node,LeafOffDiag}

# #Define partition(spine) tree
mutable struct Leaf_Spine{T <: Integer}
    m::T
    col_idx::T
    depth::T
end

mutable struct Node_Spine{T <: Integer}
    treeL::Union{Node_Spine{T},Leaf_Spine{T}}
    treeR::Union{Node_Spine{T},Leaf_Spine{T}}
    m::T
    col_idx::T
    depth::T
end

Tree = Union{Leaf_Spine,Node_Spine}

mutable struct Leaf_Gtree{T1 <: Integer, T2 <: Number}
    m::T1
    col_idx::T1
    depth::T1    
    g::Array{T2,2}
end

mutable struct Node_Gtree{T1 <: Integer, T2 <: Number}
    gtreeL::Union{Node_Gtree,Leaf_Gtree}
    gtreeR::Union{Node_Gtree,Leaf_Gtree}
    m::T1
    col_idx::T1
    depth::T1
    g::Array{T2,2}
end

Gtree = Union{Leaf_Gtree, Node_Gtree}

end

################### Shiv's Way of Defining a Generic Tree #######################
#Define g_Tree also?  Then feed g_Leaf, g_Node to Leaf_Tree and Node_Tree when I call it.

# type Leaf_Tree{T1 <: Integer, T2 <: Number, T3 <: Any}
#     depth::T1    
#     data::T3
# end

# type Node_Tree{T1 <: Integer, T2 <: Number, T3 <: Any}
#     treeL::Union{Node_tree,Leaf_tree}
#     treeR::Union{Node_tree,Leaf_tree}
#     depth::T1
#     data::T3
# end

# type g_Leaf{T1 <: Number}
#     m::T1
#     g::Array{T1,2}
# end

# type g_Node{T1 <: Number}
#     m::T1
#     col_idx::T1
#     g::Array{T1,2}
# end

