function [binaryImage ] = TextToBinaryImage( textInput )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
asc= double(textInput);
len= length(asc)
binum = dec2bin(asc,8);
sz=size(binum);
binImage= zeros(len,8);
for i=1:len
    for j=1:8
      binImage(i,j)= bin2dec(binum(i,j));
    end
end
binaryImage=binImage;
imwrite(binImage,'binaryimage.bmp');
BinaryImageToText(binaryImage)
imshow(binaryImage)


end

