function Seq = utilGenMseq(Bits)
% UTILGENMSEQ generates a binary M-sequence
% -----------------------------------------------------------------------------
% utilGenMseq 
% ===========
% Generates a Maximum Length Sequence of n bits by utilizing a linear feedback
% shift register with an XOR gate on the tap bits
%
% This function is a slightly optimised version of that written by Christopher
% Brown, cbrown@phi.luc.edu and distributed via the MATHWORKS file exchange as
% function 'mls'.
%
% Reference:
% Davies, W.D.T. (June, July, August, 1966). Generation and properties of
% maximum-length sequences. Control, 302-4, 364-5,431-3.
%
% Note: This function can take a very long time to run for longer
%       sequence lengths. For generating M-sequence visual stimuli, Setting
%       'Bits' to values of between 10 and 14 are generally appropriate for
%       average frame rates and trial lengths.
%
% Usage
% =====
% Seq = utilGenMseq(Bits);
%
% Parameters
% ==========
% Bits - The bit length of the sequence. I.e. the returned sequence length
%        would be (2^Bits)-1; This function can accept bit lengths of
%        between 2 and 24.
%
% Return Values
% =============
% Seq - A vector of 1's & 0's that is (2^n)-1 in length.
%
% Reference page in Help browser
% <a href="matlab:web(['jar:file:',which('crsHelp\help.jar'),'!/crs\tools\utilities\utilGenMseq.html'],'-helpbrowser');">utilGenMseq HTML help.</a>                                                                                                                                                        
%
% -----------------------------------------------------------------------------
  switch Bits    %[T1,T2,T3,T4]
  case 2;  Taps = [ 1, 2      ]; % Assign taps which will yeild a maximum
  case 3;  Taps = [ 1, 3      ]; % length sequence for a given bit length.
  case 4;  Taps = [ 1, 4      ]; % Reference: Vanderkooy, JAES, 42(4), 1994.
  case 5;  Taps = [ 2, 5      ];             
  case 6;  Taps = [ 1, 6      ];            
  case 7;  Taps = [ 1, 7      ];
  case 8;  Taps = [ 2, 3, 4, 8];
  case 9;  Taps = [ 4, 9      ];
  case 10; Taps = [ 3,10      ];
  case 11; Taps = [ 2,11      ];
  case 12; Taps = [ 1, 4, 6,12];
  case 13; Taps = [ 1, 3, 4,13];
  case 14; Taps = [ 1, 3, 5,14];
  case 15; Taps = [ 1,15      ];
  case 16; Taps = [ 2, 3, 5,16];
  case 17; Taps = [ 3,17      ];
  case 18; Taps = [ 7,18      ];
  case 19; Taps = [ 1, 2, 5,19];
  case 20; Taps = [ 3,20      ];
  case 21; Taps = [ 2,21      ];
  case 22; Taps = [ 1,22      ];
  case 23; Taps = [ 5,23      ];
  case 24; Taps = [ 1, 3, 4,24];
  otherwise; disp('\n Input Bits must be between 2 and 24'); return; end;

  rand('state',sum(100*clock));
  while(true); Buff = round(rand(1,Bits)); if find(Buff==1); break; end; end;

  SeqLen  = (2^Bits)-1;
  Seq     = ones(SeqLen,1);
  NumTaps = numel(Taps);

  if NumTaps==2
    T1 = Taps(1); T2 = Taps(2);
    for(ii=SeqLen:-1:1);
      Buff    = [xor(Buff(T1),Buff(T2)),Buff(1:Bits-1)];
      Seq(ii) = Buff(1);
    end
  elseif NumTaps==4
    T1 = Taps(1); T2 = Taps(2); T3 = Taps(3); T4 = Taps(4);
    for(ii=SeqLen:-1:1);
      Buff    = [xor(xor(Buff(T1),Buff(T2)),xor(Buff(T3),Buff(T4))),Buff(1:Bits-1)];
      Seq(ii) = Buff(1);
    end
  end