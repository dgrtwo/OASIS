def zeros((dim1, dim2)):
    """return a 2D array of zeroes with the given dimensions"""
    return [[0 for i in range(dim2)] for j in range(dim1)]

def myalign(seq1, seq2):
    """return a local alignment"""
    m,n =  len(seq1),len(seq2) #length of two sequences
    
    penalty=-4;   #define the gap penalty
    match_score_num = 1
    mismatch_score_num = -3
    
    #generate DP table and traceback path pointer matrix
    score=zeros((m+1,n+1))   #the DP table
    pointer=zeros((m+1,n+1))  #to store the traceback path
    
    P=0;
    
    def match_score(alpha,beta):
        if alpha == beta: return match_score_num
        return mismatch_score_num
    
    max_score=P;  #initial maximum score in DP table
    
    max_i = 0
    max_j = 0
    
    #calculate DP table and mark pointers
    for i in range(1,m+1):
        for j in range(1,n+1):
            score_up=score[i-1][j]+penalty;
            score_down=score[i][j-1]+penalty;
            score_diagonal=score[i-1][j-1]+match_score(seq1[i-1],seq2[j-1]);
            score[i][j]=max(0,score_up,score_down,score_diagonal);
            if score[i][j]==0:
                pointer[i][j]=0; #0 means end of the path
            if score[i][j]==score_up:
                pointer[i][j]=1; #1 means trace up
            if score[i][j]==score_down:
                pointer[i][j]=2; #2 means trace left
            if score[i][j]==score_diagonal:
                pointer[i][j]=3; #3 means trace diagonal
            
            if score[i][j]>=max_score:
                max_i=i;
                max_j=j;
                max_score=score[i][j];
    
    align1,align2='',''; #initial sequences
    
    i,j=max_i,max_j; #indices of path starting point
    
    #traceback, follow pointers
    while pointer[i][j]!=0:
        if pointer[i][j]==3:
            align1=align1+seq1[i-1];
            align2=align2+seq2[j-1];
            i=i-1;
            j=j-1;
        elif pointer[i][j]==2:
            align1=align1+'-';
            align2=align2+seq2[j-1];
            j=j-1;
        elif pointer[i][j]==1:
            align1=align1+seq1[i-1];
            align2=align2+'-';
            i=i-1;
        
    align1=align1[::-1]; #reverse sequence 1
    align2=align2[::-1]; #reverse sequence 2
    
    i,j=0,0;
    
    #calcuate identity, similarity, score and aligned sequeces
    def result(align1,align2):
        symbol='';
        found=0;
        score=0;
        penalty=-4;
        identity,similarity=0,0;
        for i in range(0,len(align1)):
            if align1[i]==align2[i]:     #if two AAs are the same, then output the letter
                symbol=symbol+align1[i];
                identity=identity+1;
                similarity=similarity+1;
                score=score+match_score(align1[i],align2[i]);
            elif align1[i]!=align2[i] and align1[i]!='-' and align2[i]!='-': #if they are not identical and none of them is gap
    
                score=score+match_score(align1[i],align2[i]);
                if found==1:
                    symbol=symbol+':';   #if they are similar AA, output ':'
                    similarity=similarity+1;
                if found==0:
                        symbol=symbol+' ';  #o/w, output a space
        
                found=0;
            elif align1[i]=='-' or align2[i]=='-':   #if one of them is a gap, output a space
                symbol=symbol+' ';
                score=score+penalty;
        
        if len(align1) == 0 or len(align2) == 0: return False, False, 0, 0, 0 
        
        identity=float(identity)/len(align1)*100;
        similarity=float(similarity)/len(align2)*100;
        
        return (align1, align2, score, identity, similarity)
        
    return result(align1,align2)