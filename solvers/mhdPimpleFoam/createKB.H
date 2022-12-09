Info<< "Creating field magnetic energy KB\n" << endl;

volScalarField KB("KB", DBU*magSqr(B));

if (B.nOldTimes())
{
    volVectorField* Bold  = &B.oldTime();
    volScalarField* KBold = &KB.oldTime();
    
    *KBold == DBU*magSqr(*Bold);
    
    while (Bold->nOldTimes())
    {
        Bold  = &Bold->oldTime();
        KBold = &KBold->oldTime();
        
        *KBold == DBU*magSqr(*Bold);
    }
}
