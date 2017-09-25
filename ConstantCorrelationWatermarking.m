% Copyright, Eze P.U., Udaya P. 2017.
% Email: peze@student.unimelb.edu.au
clear;
clc;
blocksize=8;
desiredcorrelation= 8.2521;% You can choose any value to experiment with

% [-1, 1, 1, -1, 1, -1, 1, 1
goldcode=[Your 8x8 goldcode goes here];





dicomlist = dir(fullfile(pwd,'DataSets/series2/','*.dcm'));
for cnt = 1 : numel(dicomlist)
    I{cnt} = dicomread(fullfile(pwd,'DataSets/series2/',dicomlist(cnt).name));  
end
% return;
% =======DEFINE GLOBAL VARIABLES FO RUNNING THE SIMULATION=========
withoutwatermark=[];
zeroAlpha = [];
oneAlpha=[];
extractzerocorrelation=[];% correlations to extract watermark in blocks
onesextracted=[];extractedwatermark=[];
ssm0=[];psn0=[];ssm1=[];psn1=[];
BER=0;berOverall=[];fnegative=[];bervalue = [];
onescount=0;zeroscount=0;
countzeromessage=0;countonemessage=0;
overallTotalextracted=0;grandTotalEmbed=0;
epsilon = 0.5;% Tolerance for overflow, underflow and unintentional attack
for cnt = 1 : numel(dicomlist)
    dicom=strcat('DataSets/series2/',dicomlist(cnt).name);
    % cell2mat()
     info=dicominfo(dicom);% extract the metadata from the image
     metadata = info;
     im=dicomread(info);
    %  im = I(1);
     pixeldepth=info.BitDepth;
    peakval=2^(pixeldepth);
     if length(size(im)) >2
        im=rgb2gray(im);
     end
    [row,col]=size(im);
    im_uncast = im;%NOT cast from int16 to uint16
     im = uint16(im); % Necessary to ensure pixel values are not accepted at underflow and overflow values
    

    block=im(1:8,1:8);
    [N, M]=size(block);
    message=round(rand(row/8,col/8));% Experimental message from random bits

   
    

    %==========ESTIMATION OF LINEAR CORRELATION THRESHOLD AND RECOMMENDED ALPHA
    % ========= CORRELATION CONSTANT FOR EACH BLOCK



    i=1;j=1;
    caption=strcat('Watermark Embedding and Extraction for:', dicom,'. Please wait...');
    display(caption);
    bitToEmbed=1;
    finalstego = zeros(row,col);
    modifiedcover = zeros(row,col);
    totalmodifiedpixels=0;
    none=0;
    x=1;y=1;
    errorbars =[];
    ber=0; % error bits
    falsenegative =0;
    while i<row
        while j<col
            block=im(i:i+7,j:j+7);
            withoutwatermark=[withoutwatermark linearCorrelation(block,goldcode)];
            bitToEmbed=message(x,y);

            alpha=computeAlpha(block,desiredcorrelation,bitToEmbed,goldcode);

    %             alpha=desiredcorrelation;
            %==== HISTOGRAM EQUALIZATION TO OVERCOME OVERFLOW AND UNDERFLOW
            [modpixels,block]=HistogramEqualize(block,pixeldepth,alpha);%Equalise histogram to remove underflow and overflow
            modifiedcover(i:i+7,j:j+7)=block;
            totalmodifiedpixels=totalmodifiedpixels + modpixels;
    %         alpha=computeAlpha(block,desiredcorrelation,bitToEmbed,goldcode);
            %================Ends Here===================%
            zeroAlpha=[zeroAlpha alpha];
            watermark=(alpha.*goldcode);
            if (bitToEmbed==0) % Embed 0-bit according to our paper
                stego=(double(block)+ watermark);
                countzeromessage=countzeromessage+1;
            else
                stego=(double(block)- watermark);% Embed 1-bit according to our paper
               countonemessage=countonemessage+1;
            end

            SSM = ssim(uint16(stego),block);% more than 8-bit depth
            PSN = psnr(uint16(stego),block,peakval);
    %         PSN = psnr(uint8(stego),block);
    %         SSM = ssim(uint8(stego),block);
            ssm0=[ssm0 SSM];psn0=[psn0 PSN];
    %         normalisedstego = (uint16(stego) - mean2(uint16(stego)))/std2(uint16(stego))
            corr0=linearCorrelation(uint16(stego),goldcode);
    %         corr0=linearCorrelation(uint8(stego),goldcode);
            finalstego(i:i+7,j:j+7)=stego;
            extractzerocorrelation= [extractzerocorrelation corr0];
            errorbars =[errorbars abs(corr0)-desiredcorrelation];
    

            % Modified criteria for extraction and tamper detection
            if corr0 < (-desiredcorrelation + epsilon) && (corr0 > -desiredcorrelation - epsilon) % corr = -p+-0.5
                extractedwatermark=[extractedwatermark bitToEmbed];
                onescount=onescount+1; 
                 if message(x,y)~= 1
                   ber = ber+1;
                 end
            elseif corr0 < (desiredcorrelation + epsilon) && (corr0 > desiredcorrelation - epsilon) % corr = p+-0.5
                extractedwatermark=[extractedwatermark bitToEmbed];
                 zeroscount=zeroscount+1; 
               if message(x,y)~= 0
                   ber = ber+1;
               end
            else
                extractedwatermark=[extractedwatermark 9];%False negative value - non-binary
                falsenegative =falsenegative+1;
                ber = ber+1;
                finalstego(i:i+7,j:j+7)=2^(pixeldepth)-1; % Flagging potentially tampered blocks
            end


            % Increment to next bit in message matrix
            y=y+1;
            if y>(col/N)
                x=x+1;
                y=1;
            end
            if x>(row/M)
                x=1;
            end
   
            j=j+8;%Go to next column
        end
        j=1;
        i=i+8;% Go to next row
    end
    totalembed=countonemessage+countzeromessage;
    countonemessage=0;countzeromessage=0; % Reset counters for next image
    grandTotalEmbed = grandTotalEmbed +totalembed;
    totalextracted=zeroscount+onescount;
    zeroscount=0;onescount=0;
    overallTotalextracted = overallTotalextracted + totalextracted;
    totlalextracted=0;
    bervalue=[bervalue ber];
    BER = (ber)/(totalembed);
    ber=0;
    falsenegative=0;
    [r,c]=size(message);
    
    berOverall=[berOverall BER];fnegative=[fnegative falsenegative];
    BER=0;
    filename=strcat('watermarkedDicom',num2str(cnt),'.dcm')
    falseNegative = (falsenegative)/(r*c);
    %dicomwrite(X, 'ct_copy.dcm', metadata, 'CreateMode', 'copy');
    dicomwrite(uint16(finalstego),filename, metadata);

end % End for EACH image

ber_average = sum(bervalue)/grandTotalEmbed
bervalue
berOverall
% Overall unwatermarked statistics
 meanTg = mean(withoutwatermark)% Mean Correlation without embedded data
%
 variancethreshold=var(withoutwatermark);
 standarddev= std(withoutwatermark);
 recommendedConstantCorrelation=meanTg+standarddev; % use this as p



figure
plot(zeroAlpha,psn0,'*')
xlabel('Embedding Strength,a')
ylabel('PSNR (dB)')
caption=strcat('Embedding Strength, a Vs PSNR at p=',num2str(desiredcorrelation));
title( caption)
grid on
grid minor

figure
plot(zeroAlpha,ssm0,'o')
xlabel('Embedding Strength,a')
ylabel('SSIM')
caption=strcat('Embedding Strength, a Vs SSIM at p=',num2str(desiredcorrelation));
title(caption)
grid on
grid minor

figure
bins=30;
histogram(withoutwatermark,bins)
% axis([-bins/2 bins/2 0 1500])
title('Correlation values for Unwatermarked Image blocks')
xlabel('Linear Correlation')
ylabel('Number of Image Blocks')
grid on
grid minor



figure
bins=30;
histogram(extractzerocorrelation,bins)

caption=strcat('Deviations from constant correlation, p=',num2str(desiredcorrelation));
title(caption)
xlabel('Linear Correlation')
ylabel('Number of Image Blocks')
grid on
grid minor

figure
bins=30;
histogram(zeroAlpha,bins)

caption=strcat('Distribution of Embedding Strength, a at p=',num2str(desiredcorrelation));
title(caption)
xlabel('Embedding strength,a')
ylabel('Number of Image Blocks')
grid on
grid minor

figure
bins=30;
histogram(ssm0,bins)
caption=strcat('Distribution of SSIM at p=',num2str(desiredcorrelation));
title(caption)
xlabel('SSIM')
ylabel('Number of Image Blocks')
grid on
grid minor

figure
bins=30;
histogram(psn0,bins)
caption=strcat('Distribution of PSNR at p=',num2str(desiredcorrelation));
title(caption)
xlabel('PSNR(dB)')
ylabel('Number of Image Blocks')
grid on
grid minor

figure
bins=20;
histogram(berOverall)
caption=strcat('Distribution of BER at p=',num2str(desiredcorrelation));
title(caption)
xlabel('BER')
ylabel('Number of Images')
grid on
grid minor

figure
bins=20;
histogram(bervalue)
caption=strcat('Distribution of BER at p=',num2str(desiredcorrelation));
title(caption)
xlabel('BER (No. of Bits)')
ylabel('Number of Images')
grid on
grid minor


figure
subplot(2,2,1)
% imshow(im)
imshow(dicom,'DisplayRange',[])
title('Original Cover')

subplot(2,2,2)
imshow(uint16(finalstego),'DisplayRange',[])
% imshow(coverstego)
caption=strcat('Stego Image at p=',num2str(desiredcorrelation));
title(caption)

subplot(2,2,3)
imhist(im)
title('Histogram of Cover Image')

subplot(2,2,4)
% imhist(uint8(finalstego))
imhist(uint16(finalstego))
caption=strcat('Histogram of Stego Image at p=',num2str(desiredcorrelation));
title(caption)




