function [ textOutput ] = BinaryImageToText( binaryImage )
% BinaryImageToText - Converts Nx8 binary image to ASCII code Texts
% binaryImage - Nx8 array of binary image
% textOutput - The equivalent characters of the binary image
% Convert binary image to text
binstring ='';
sz= size(binaryImage);
N = sz(1);
len = N;
binarystring = zeros(1,len);
for i=1:len
    for j=1:8
        if binaryImage(i,j)==1
            binstring = strcat(binstring,'1');
        else
            binstring = strcat(binstring,'0'); 
        end
    end
    binarystring(1,i) = bin2dec(binstring);
    binstring='';
end
textOutput=char(binarystring);

end

