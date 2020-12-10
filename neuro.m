function R=neuro(V,S,Inp)

R=max(0,Inp'*reshape(V(1:S(1)*S(2)),S(1),S(2)))*reshape(V(S(1)*S(2)+1:end),S(2),S(3));

end