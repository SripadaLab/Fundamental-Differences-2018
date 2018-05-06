function [Q] = calc_modularity_und(A,C,method,gamma)
%CALC_MODULARITY_UND    Calculate modularity for a fixed community structure
%
% This function is adapted from "modularity_und" available from the
% Brain Connectivity Toolbox
% (https://sites.google.com/site/bctnet/). Whereas modularity_und
% learns a community structure that optimizes the modularity Q,
% this function simply evaluates the modularity for a given
% community structure C.

if ~exist('gamma','var')
    gamma = 1;
end

if ~exist('method','var')
    method = 'bct';
end

switch method
  case 'bct'

    N=length(A);                            %number of vertices
    K=sum(A);                               %degree
    m=sum(K);                               %number of edges (each undirected edge is counted twice)
    B=A-gamma*(K.'*K)/m;                    %modularity matrix



    s=C(:,ones(1,N));                      %compute modularity
    Q=~(s-s.').*B/m;                       % intracomm mask times B/m
    Q=sum(Q(:));
    
  case 'gomez'
    
    wij_p = max(0,A);
    wij_n = max(0,-A);

    wi_p = sum(wij_p,1);
    wi_n = sum(wij_n,1);

    wp2 = sum(wi_p);
    wn2 = sum(wi_n);

    q = 0;
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            if C(i)==C(j)
                q = q + A(i,j) - (wi_p(i) * wi_p(j) / wp2 - wi_n(i) * wi_n(j) / wn2);
            end
        end
    end

    Q = q / (wp2 + wn2);

    case 'newman'
        ki = sum(A,1);
        l = sum(ki);

        q = 0;
        for i = 1:size(A,1)
            for j = 1:size(A,1)
                if C(i)==C(j)
                    q = q  + A(i,j) - ki(i)* ki(j) / l;
                end
            end
        end

        Q = q / l;
end
