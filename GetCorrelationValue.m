function corr_val=GetCorrelationValue(x,y)

    corr_val=dot(x, y)/(norm(x,2)*norm(y,2));

end