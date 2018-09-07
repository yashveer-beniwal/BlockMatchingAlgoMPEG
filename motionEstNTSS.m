% Computes motion vectors using *NEW* Three Step Search method
%
% Based on the paper by R. Li, b. Zeng, and M. L. Liou
% IEEE Trans. on Circuits and Systems for Video Technology
% Volume 4, Number 4, August 1994 :  Pages 438:442
%
% Input
%   imgP : The image for which we want to find motion vectors
%   imgI : The reference image
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   motionVect : the motion vectors for each integral macroblock in imgP
%   NTSScomputations: The average number of points searched for a macroblock
%
% Written by Aroh Barjatya


function [motionVect, NTSScomputations] = motionEstNTSS(imgP, imgI, mbSize, p)

[row col] = size(imgI);

vectors = zeros(2,row*col/mbSize^2);
costs = ones(3, 3) * 65537;


% we now take effectively log to the base 2 of p
% this will give us the number of steps required

L = floor(log10(p+1)/log10(2));   
stepMax = 2^(L-1);

computations = 0;

% we start off from the top left of the image
% we will walk in steps of mbSize
% for every marcoblock that we look at we will look for
% a close match p pixels on the left, right, top and bottom of it

mbCount = 1;
for i = 1 : mbSize : row-mbSize+1
    for j = 1 : mbSize : col-mbSize+1
        
        % the NEW three step search starts


        
        x = j;
        y = i;
        
        % In order to avoid calculating the center point of the search
        % again and again we always store the value for it from the
        % previous run. For the first iteration we store this value outside
        % the for loop, but for subsequent iterations we store the cost at
        % the point where we are going to shift our root.
        %
        % For the NTSS, we find the minimum first in the far away points
        % we then find the minimum for the close up points
        % we then compare the minimums and which ever is the lowest is where
        % we shift our root of search. If the minimum is the center of the
        % current window then we stop the search. If its one of the
        % immediate close to the center then we will do the second step
        % stop. And if its in the far away points, then we go doing about
        % the normal TSS approach
        % 
        % more details in the code below or read the paper/literature
        
        costs(2,2) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                    imgI(i:i+mbSize-1,j:j+mbSize-1),mbSize);
        stepSize = stepMax; 
        computations = computations + 1;

        % This is the calculation of the outer 8 points
        % m is row(vertical) index
        % n is col(horizontal) index
        % this means we are scanning in raster order
        for m = -stepSize : stepSize : stepSize        
            for n = -stepSize : stepSize : stepSize
                refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                refBlkHor = x + n;   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                     || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                     continue;
                end

                costRow = m/stepSize + 2;
                costCol = n/stepSize + 2;
                if (costRow == 2 && costCol == 2)
                    continue
                end
                costs(costRow, costCol ) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                    imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                computations = computations + 1;
            end
        end
        
        % Now we find the vector where the cost is minimum
        % and store it ... 
        
        [dx, dy, min1] = minCost(costs);      % finds which macroblock in imgI gave us min Cost
            
              
        % Find the exact co-ordinates of this point

        x1 = x + (dx-2)*stepSize;
        y1 = y + (dy-2)*stepSize;
            
        % Now find the costs at 8 points right next to the center point
        % (x,y) still points to the center
        
        stepSize = 1;
        for m = -stepSize : stepSize : stepSize        
            for n = -stepSize : stepSize : stepSize
                refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                refBlkHor = x + n;   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                     || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                     continue;
                end

                costRow = m/stepSize + 2;
                costCol = n/stepSize + 2;
                if (costRow == 2 && costCol == 2)
                    continue
                end
                costs(costRow, costCol ) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                    imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                computations = computations + 1;
            end
        end
        
        % now find the minimum amongst this
        
        [dx, dy, min2] = minCost(costs);      % finds which macroblock in imgI gave us min Cost
            
              
        % Find the exact co-ordinates of this point

        x2 = x + (dx-2)*stepSize;
        y2 = y + (dy-2)*stepSize;
        
        % the only place x1 == x2 and y1 == y2 will take place will be the
        % center of the search region
        
        if (x1 == x2 && y1 == y2)
            % then x and y still remain pointing to j and i;
            NTSSFlag = -1; % this flag will take us out of any more computations 
        elseif (min2 <= min1)
            x = x2;
            y = y2;
            NTSSFlag = 1; % this flag signifies we are going to go into NTSS mode
        else
            x = x1;
            y = y1;
            NTSSFlag = 0; % This value of flag says, we go into normal TSS
        end
        
        
        if (NTSSFlag == 1)
            % Now in order to make sure that we dont calcylate the same
            % points again which were in the initial center window we take
            % care as follows
            
            costs = ones(3,3) * 65537;
            costs(2,2) = min2;
            stepSize = 1;
            for m = -stepSize : stepSize : stepSize        
                for n = -stepSize : stepSize : stepSize
                    refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                    refBlkHor = x + n;   % col/Horizontal co-ordinate
                    if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                           || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                        continue;
                    end
                    
                    if ( (refBlkVer >= i - 1  && refBlkVer <= i + 1) ...
                            && (refBlkHor >= j - 1  && refBlkHor <= j + 1) )
                        continue;
                    end
                    
                    costRow = m/stepSize + 2;
                    costCol = n/stepSize + 2;
                    if (costRow == 2 && costCol == 2)
                        continue
                    end
                    costs(costRow, costCol ) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                         imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                    computations = computations + 1;
                end
            end
                
            % now find the minimum amongst this
        
            [dx, dy, min2] = minCost(costs);      % finds which macroblock in imgI gave us min Cost
            
            % Find the exact co-ordinates of this point and stop

            x = x + (dx-2)*stepSize;
            y = y + (dy-2)*stepSize;            
            
        elseif (NTSSFlag == 0)
            % this is when we are going about doing normal TSS business
            costs = ones(3,3) * 65537;
            costs(2,2) = min1;
            stepSize = stepMax / 2;
            while(stepSize >= 1)  
                for m = -stepSize : stepSize : stepSize        
                    for n = -stepSize : stepSize : stepSize
                        refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                        refBlkHor = x + n;   % col/Horizontal co-ordinate
                        if ( refBlkVer < 1 || refBlkVer+mbSize-1 > row ...
                            || refBlkHor < 1 || refBlkHor+mbSize-1 > col)
                            continue;
                        end

                        costRow = m/stepSize + 2;
                        costCol = n/stepSize + 2;
                        if (costRow == 2 && costCol == 2)
                            continue
                        end
                        costs(costRow, costCol ) = costFuncMAD(imgP(i:i+mbSize-1,j:j+mbSize-1), ...
                                imgI(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                        computations = computations + 1;
                    
                    end
                end
        
                % Now we find the vector where the cost is minimum
                % and store it ... this is what will be passed back.
        
                [dx, dy, min] = minCost(costs);      % finds which macroblock in imgI gave us min Cost
            
            
                % shift the root for search window to new minima point

                x = x + (dx-2)*stepSize;
                y = y + (dy-2)*stepSize;
            
                stepSize = stepSize / 2;
                costs(2,2) = costs(dy,dx);
            
            end
        end

        vectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        vectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        mbCount = mbCount + 1;
        costs = ones(3,3) * 65537;
    end
end

motionVect = vectors;
NTSScomputations = computations/(mbCount - 1);
                    