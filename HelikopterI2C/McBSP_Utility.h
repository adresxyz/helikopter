#ifndef _McBSP_Utility_H
#define _McBSP_Utility_H

/*McBSP Handlers*/
static MCBSP_Handle hMcbsp1;  		/*Handler for MCBSP_DEV1*/
static MCBSP_Handle hMcBSP_ADC; 	/*Handler for MCBSP_DEV0 ADC Communication*/
static MCBSP_Handle hMcBSP_DAC; 	/*Handler for MCBSP_DEV1 DAC Communication*/

/*McBSP IRQ Events IDs */
static Uint32 Xmt1_EventId; 		/*hMcbsp1 Transmit Event ID*/
static Uint32 Recv1_EventId; 		/*hMcbsp1 Receive Event ID*/

/*Help Functions*/
void Polling_Transmit(MCBSP_Handle _hMcbsp, int _WriteValue);
void Polling_Single_Receive(MCBSP_Handle _hMcbsp, int* _ReadTo);
void MCBSP_Close_All();

/*McBSP Configuration Structures */

/*--------------------------------------------------McBSP ADC Confing----------------------------------------------*/


static MCBSP_Config McBSP_ADC_Confing = {
MCBSP_SPCR_RMK( /* Serial Port Control Reg. (SPCR) */

		MCBSP_SPCR_FREE_YES, /* Serial clock free running mode(FREE)
		                              MCBSP_SPCR_FREE_NO - During emulation halt,SOFT bit
		                                   determines operation of McBSP.
		                              MCBSP_SPCR_FREE_YES - During emulation halt, serial
		                                   clocks continue to run.                     */

		MCBSP_SPCR_SOFT_YES, /* Serial clock emulation mode(SOFT)
		                              MCBSP_SPCR_SOFT_NO -  In conjunction with FREE = 0,
		                                   serial port clock stops immediately during
		                                   emulation halt,thus aborting any transmissions.
		                              MCBSP_SPCR_SOFT_YES - In conjunction with FREE = 0,
		                                   during emulation halt, serial port clock stops
		                                   after completion of current transmission.   */

		MCBSP_SPCR_FRST_YES, /* Frame sync generator reset(FRST)
		                              MCBSP_SPCR_FRST_YES /MCBSP_SPCR_FRST_0 - Frame sync
		                                   generation logic is reset.
		                              MCBSP_SPCR_FRST_NO /MCBSP_SPCR_FRST_1  - Frame sync
		                                   signal is generated after eight CLKG clocks. */

		MCBSP_SPCR_GRST_NO, /* Sample rate generator reset(GRST)
		                             MCBSP_SPCR_GRST_YES /MCBSP_SPCR_GRST_0  - Reset
		                             MCBSP_SPCR_GRST_NO /MCBSP_SPCR_GRST_1 -Out of reset*/

		MCBSP_SPCR_XINTM_XRDY,/* Transmit interrupt mode(XINTM)
		                             MCBSP_SPCR_XINTM_XRDY - XINT driven by XRDY
		                             MCBSP_SPCR_XINTM_EOS  - XINT generated by
		                                   end-of-subframe in multichannel operation
		                             MCBSP_SPCR_XINTM_FRM  - XINT generated by a new frame
		                                   synchronization
		                             MCBSP_SPCR_XINTM_XSYNCERR - XINT generated by XSYNCERR
		                                                                                */
		MCBSP_SPCR_XSYNCERR_NO,/* Transmit synchronization error
		                             MCBSP_SPCR_XSYNCERR_NO /MCBSP_SPCR_XSYNCERR_0  - No
		                                   frame synchronization error
		                             MCBSP_SPCR_XSYNCERR_YES /MCBSP_SPCR_XSYNCERR_1 - Frame
		                                   synchronization error detected by McBSP       */

		MCBSP_SPCR_XRST_YES, /* Transmitter reset(XRST). This resets or enables
		                            transmitter.
		                             MCBSP_SPCR_XRST_YES /MCBSP_SPCR_XRST_0  - Serial port
		                                   transmitter is disabled and is in reset state.
		                             MCBSP_SPCR_XRST_NO /MCBSP_SPCR_XRST_1   - Serial port
		                                   transmitter is enabled.                       */

		MCBSP_SPCR_DLB_OFF,   /* Digital loopback(DLB) mode
		                             MCBSP_SPCR_DLB_OFF   - DLB Disabled
		                             MCBSP_SPCR_DLB_ON    - DLB Enabled                  */

		MCBSP_SPCR_RJUST_RZF,/*Receive data sign-extension and justification mode(RJUST)
		                             MCBSP_SPCR_RJUST_RZF - Right-justify and zero-fill
		                                   MSBs in DRR.
		                             MCBSP_SPCR_RJUST_RSE - Right-justify and sign-extend
		                                   MSBs in DRR.
		                             MCBSP_SPCR_RJUST_LZF - Left-justify and zero-fill
		                                   LSBs in DRR.                                  */

		MCBSP_SPCR_CLKSTP_DISABLE,/* Clock stop(CLKSTP) mode
		                             MCBSP_SPCR_CLKSTP_DISABLE - Disabled. Normal
		                                  clocking enabled for non-SPI mode.
		                             MCBSP_SPCR_CLKSTP_NODELAY -Clock starts without delay.
		                             MCBSP_SPCR_CLKSTP_DELAY   - Clock starts with delay.*/

		MCBSP_SPCR_DXENA_OFF,/* DX Enabler(DXENA) -Extra delay for DX turn-on time.
		                             Only first bit of data is delayed.
		                             MCBSP_SPCR_DXENA_OFF - DX enabler is off.
		                             MCBSP_SPCR_DXENA_ON  - DX enabler is on.            */

		MCBSP_SPCR_RINTM_RRDY,/* Receive interrupt(RINT) mode
		                             MCBSP_SPCR_RINTM_RRDY - RINT driven by RRDY
		                             MCBSP_SPCR_RINTM_EOS  - RINT generated by
		                                  end-of-subframe in multichannel operation
		                             MCBSP_SPCR_RINTM_FRM  - RINT generated by a new frame
		                                  synchronization
		                             MCBSP_SPCR_RINTM_RSYNCERR -RINT generated by RSYNCERR*/

		MCBSP_SPCR_RSYNCERR_NO,/* Receive synchronization error(RSYNCERR)
		                             MCBSP_SPCR_RSYNCERR_NO /MCBSP_SPCR_RSYNCERR_0 - No
		                                  frame synchronization error
		                             MCBSP_SPCR_RSYNCERR_YES /MCBSP_SPCR_RSYNCERR_1 - Frame
		                                  synchronization error detected by McBSP         */

		MCBSP_SPCR_RRST_NO /* Receiver reset(RRST). This resets or enables
		                           receiver.
		                             MCBSP_SPCR_RRST_YES /MCBSP_SPCR_RRST_0 - Serial
		                                 port receiver is disabled and is in reset state.
		                             MCBSP_SPCR_RRST_NO  /MCBSP_SPCR_RRST_1 - Serial
		                                  port receiver is enabled.                    */

),
MCBSP_RCR_RMK(  /* Receiver Control Reg. (RCR) */

		MCBSP_RCR_RPHASE_SINGLE, /* Receive phases
                             MCBSP_RCR_RPHASE_SINGLE - Single phase frame
                             MCBSP_RCR_RPHASE_DUAL   - Dual phase frame        */

		MCBSP_RCR_RFRLEN2_OF(0), /* Receive frame length in phase 2(RFRLEN2)
                              000 0000b: 1 word per phase
                              000 0001b: 2 words per phase
                              . . . . . . . . . . . .
                              111 1111b: 128 words per phase                   */

		MCBSP_RCR_RWDLEN2_20BIT,/* Receive element length in phase 2(RWDLEN2)
                             MCBSP_RCR_RWDLEN2_8BIT  -  8  bits
                             MCBSP_RCR_RWDLEN2_12BIT -  12 bits
                             MCBSP_RCR_RWDLEN2_16BIT -  16 bits
                             MCBSP_RCR_RWDLEN2_20BIT -  20 bits
                             MCBSP_RCR_RWDLEN2_24BIT -  24 bits
                             MCBSP_RCR_RWDLEN2_32BIT -  32 bits                */

		MCBSP_RCR_RCOMPAND_MSB,/* Receive companding mode (RCOMPAND)
                             MCBSP_RCR_RCOMPAND_MSB  - No companding.Data
                                  transfer starts with MSB first.
                             MCBSP_RCR_RCOMPAND_8BITLSB - No companding,
                                  8-bit data. Transfer starts with LSB first.
                                  Applicable to 8-bit data or 32-bit data in
                                  data reversal mode.
                             MCBSP_RCR_RCOMPAND_ULAW -  Compand using m-law for
                                  receive data. Applicable to 8-bit data only
                             MCBSP_RCR_RCOMPAND_ALAW - Compand using A-law for
                                  receive data. Applicable to 8-bit data only  */

		MCBSP_RCR_RFIG_YES, /* Receive frame ignore(RFIG)
                             MCBSP_RCR_RFIG_NO  - Unexpected receive frame
                                  synchronization pulses restart the transfer.
                             MCBSP_RCR_RFIG_YES - Unexpected receive frame
                                  synchronization pulses are ignored.          */

		MCBSP_RCR_RDATDLY_1BIT,/* Receive data delay(RDATDLY)
                             MCBSP_RCR_RDATDLY_0BIT - 0 bit data delay
                             MCBSP_RCR_RDATDLY_1BIT - 1 bit data delay
                             MCBSP_RCR_RDATDLY_2BIT - 2 bit data delay         */

		MCBSP_RCR_RFRLEN1_OF(0), /* Receive frame length in phase 1(RFRLEN1)
                             000 0000b: 1 word per phase
                             000 0001b: 2 words per phase
                             . . . . . . . . . . . .
                             111 1111b: 128 words per phase                    */

		MCBSP_RCR_RWDLEN1_20BIT,/* Receive element length in phase 1(RWDLEN1)

                             MCBSP_RCR_RWDLEN1_8BIT  -  8 bits
                             MCBSP_RCR_RWDLEN1_12BIT - 12 bits
                             MCBSP_RCR_RWDLEN1_16BIT - 16 bits
                             MCBSP_RCR_RWDLEN1_20BIT - 20 bits
                             MCBSP_RCR_RWDLEN1_24BIT - 24 bits
                             MCBSP_RCR_RWDLEN1_32BIT - 32 bits                 */

		MCBSP_RCR_RWDREVRS_DISABLE/* Receive 32-bit bit reversal feature.(RWDREVRS)
                             MCBSP_RCR_RWDREVRS_DISABLE -32 bit reversal disabled
                             MCBSP_RCR_RWDREVRS_ENABLE  -32 bit reversal enabled.
                                  32-bit data is received LSB first.
                                  RWDLEN should be set for 32-bit operation.
                                  RCOMPAND should be set to 01b else operation
                                  is undefined.                                */
 ),

0x00010000, /* Transmitter Control Reg. (XCR) */

MCBSP_SRGR_RMK( /* Sample�Rate Generator Reg. (SRGR) */

		MCBSP_SRGR_GSYNC_FREE,/* Sample rate generator clock synchronization(GSYNC).
                           MCBSP_SRGR_GSYNC_FREE - The sample rate generator
                                clock CLKG) is free running.
                           MCBSP_SRGR_GSYNC_SYNC - (CLKG) is running but is
                                 resynchronized, and the frame sync signal
                                 (FSG)is generated only after the receive
                                 frame synchronization signal(FSR)is detected.
                                 Also,the frame period (FPER) is a don�t care
                                 because the period is dictated by the external
                                 frame sync pulse.                            */

		MCBSP_SRGR_CLKSP_RISING,/* CLKS polarity clock edge select(CLKSP)
                            MCBSP_SRGR_CLKSP_RISING  - The rising edge of CLKS
                                 generates CLKG and FSG.
                            MCBSP_SRGR_CLKSP_FALLING - The falling edge of CLKS
                                 generates CLKG and FSG.                      */
		MCBSP_SRGR_CLKSM_INTERNAL,/* MCBSP sample rate generator clock mode(CLKSM)
                            MCBSP_SRGR_CLKSM_CLKS   - The sample rate generator
                                 clock is derived from CLKS.
                            MCBSP_SRGR_CLKSM_INTERNAL - (Default value) The
                                 sample rate generator clock is derived from
                                 the internal clock source.                   */

		MCBSP_SRGR_FSGM_DXR2XSR,/*Sample rate generator transmit frame synchronization
                             mode.(FSGM)
                            MCBSP_SRGR_FSGM_DXR2XSR  - The transmit frame sync
                                 signal (FSX) is generated on every DXR to XSR
                                 copy.
                            MCBSP_SRGR_FSGM_FSG      - The transmit frame sync
                                 signal is driven by the sample rate generator
                                 frame sync signal, FSG.                      */

		MCBSP_SRGR_FPER_OF(21),/* Frame period(FPER)
                            Valid values: 0 to 4095                           */

		MCBSP_SRGR_FWID_OF(0),/* Frame width(FWID)
                            Valid values: 0 to 255                            */

		MCBSP_SRGR_CLKGDV_OF(5)/* Sample rate generator clock divider(CLKGDV)
                            Valid values: 0 to 255                            */

),
0x00000000, /* Multichannel Control Reg. (MCR) */
0x00000000, /* Receiver Channel Enable(RCER) */
0x00000000, /* Transmitter Channel Enable(XCER) */
0x00000500 /* Pin Control Reg. (PCR) */
};

static MCBSP_Config McBSP_DAC_Confing = {
MCBSP_SPCR_RMK( /* Serial Port Control Reg. (SPCR) */

		MCBSP_SPCR_FREE_YES, /* Serial clock free running mode(FREE)
		                              MCBSP_SPCR_FREE_NO - During emulation halt,SOFT bit
		                                   determines operation of McBSP.
		                              MCBSP_SPCR_FREE_YES - During emulation halt, serial
		                                   clocks continue to run.                     */

		MCBSP_SPCR_SOFT_YES, /* Serial clock emulation mode(SOFT)
		                              MCBSP_SPCR_SOFT_NO -  In conjunction with FREE = 0,
		                                   serial port clock stops immediately during
		                                   emulation halt,thus aborting any transmissions.
		                              MCBSP_SPCR_SOFT_YES - In conjunction with FREE = 0,
		                                   during emulation halt, serial port clock stops
		                                   after completion of current transmission.   */

		MCBSP_SPCR_FRST_YES, /* Frame sync generator reset(FRST)
		                              MCBSP_SPCR_FRST_YES /MCBSP_SPCR_FRST_0 - Frame sync
		                                   generation logic is reset.
		                              MCBSP_SPCR_FRST_NO /MCBSP_SPCR_FRST_1  - Frame sync
		                                   signal is generated after eight CLKG clocks. */

		MCBSP_SPCR_GRST_NO, /* Sample rate generator reset(GRST)
		                             MCBSP_SPCR_GRST_YES /MCBSP_SPCR_GRST_0  - Reset
		                             MCBSP_SPCR_GRST_NO /MCBSP_SPCR_GRST_1 -Out of reset*/

		MCBSP_SPCR_XINTM_XRDY,/* Transmit interrupt mode(XINTM)
		                             MCBSP_SPCR_XINTM_XRDY - XINT driven by XRDY
		                             MCBSP_SPCR_XINTM_EOS  - XINT generated by
		                                   end-of-subframe in multichannel operation
		                             MCBSP_SPCR_XINTM_FRM  - XINT generated by a new frame
		                                   synchronization
		                             MCBSP_SPCR_XINTM_XSYNCERR - XINT generated by XSYNCERR
		                                                                                */
		MCBSP_SPCR_XSYNCERR_NO,/* Transmit synchronization error
		                             MCBSP_SPCR_XSYNCERR_NO /MCBSP_SPCR_XSYNCERR_0  - No
		                                   frame synchronization error
		                             MCBSP_SPCR_XSYNCERR_YES /MCBSP_SPCR_XSYNCERR_1 - Frame
		                                   synchronization error detected by McBSP       */

		MCBSP_SPCR_XRST_NO, /* Transmitter reset(XRST). This resets or enables
		                            transmitter.
		                             MCBSP_SPCR_XRST_YES /MCBSP_SPCR_XRST_0  - Serial port
		                                   transmitter is disabled and is in reset state.
		                             MCBSP_SPCR_XRST_NO /MCBSP_SPCR_XRST_1   - Serial port
		                                   transmitter is enabled.                       */

		MCBSP_SPCR_DLB_OFF,   /* Digital loopback(DLB) mode
		                             MCBSP_SPCR_DLB_OFF   - DLB Disabled
		                             MCBSP_SPCR_DLB_ON    - DLB Enabled                  */

		MCBSP_SPCR_RJUST_RZF,/*Receive data sign-extension and justification mode(RJUST)
		                             MCBSP_SPCR_RJUST_RZF - Right-justify and zero-fill
		                                   MSBs in DRR.
		                             MCBSP_SPCR_RJUST_RSE - Right-justify and sign-extend
		                                   MSBs in DRR.
		                             MCBSP_SPCR_RJUST_LZF - Left-justify and zero-fill
		                                   LSBs in DRR.                                  */

		MCBSP_SPCR_CLKSTP_DISABLE,/* Clock stop(CLKSTP) mode
		                             MCBSP_SPCR_CLKSTP_DISABLE - Disabled. Normal
		                                  clocking enabled for non-SPI mode.
		                             MCBSP_SPCR_CLKSTP_NODELAY -Clock starts without delay.
		                             MCBSP_SPCR_CLKSTP_DELAY   - Clock starts with delay.*/

		MCBSP_SPCR_DXENA_OFF,/* DX Enabler(DXENA) -Extra delay for DX turn-on time.
		                             Only first bit of data is delayed.
		                             MCBSP_SPCR_DXENA_OFF - DX enabler is off.
		                             MCBSP_SPCR_DXENA_ON  - DX enabler is on.            */

		MCBSP_SPCR_RINTM_RRDY,/* Receive interrupt(RINT) mode
		                             MCBSP_SPCR_RINTM_RRDY - RINT driven by RRDY
		                             MCBSP_SPCR_RINTM_EOS  - RINT generated by
		                                  end-of-subframe in multichannel operation
		                             MCBSP_SPCR_RINTM_FRM  - RINT generated by a new frame
		                                  synchronization
		                             MCBSP_SPCR_RINTM_RSYNCERR -RINT generated by RSYNCERR*/

		MCBSP_SPCR_RSYNCERR_NO,/* Receive synchronization error(RSYNCERR)
		                             MCBSP_SPCR_RSYNCERR_NO /MCBSP_SPCR_RSYNCERR_0 - No
		                                  frame synchronization error
		                             MCBSP_SPCR_RSYNCERR_YES /MCBSP_SPCR_RSYNCERR_1 - Frame
		                                  synchronization error detected by McBSP         */

		MCBSP_SPCR_RRST_YES /* Receiver reset(RRST). This resets or enables
		                           receiver.
		                             MCBSP_SPCR_RRST_YES /MCBSP_SPCR_RRST_0 - Serial
		                                 port receiver is disabled and is in reset state.
		                             MCBSP_SPCR_RRST_NO  /MCBSP_SPCR_RRST_1 - Serial
		                                  port receiver is enabled.                    */

),
MCBSP_RCR_RMK(  /* Receiver Control Reg. (RCR) */

		MCBSP_RCR_RPHASE_SINGLE, /* Receive phases
                             MCBSP_RCR_RPHASE_SINGLE - Single phase frame
                             MCBSP_RCR_RPHASE_DUAL   - Dual phase frame        */

		MCBSP_RCR_RFRLEN2_OF(0), /* Receive frame length in phase 2(RFRLEN2)
                              000 0000b: 1 word per phase
                              000 0001b: 2 words per phase
                              . . . . . . . . . . . .
                              111 1111b: 128 words per phase                   */

		MCBSP_RCR_RWDLEN2_24BIT,/* Receive element length in phase 2(RWDLEN2)
                             MCBSP_RCR_RWDLEN2_8BIT  -  8  bits
                             MCBSP_RCR_RWDLEN2_12BIT -  12 bits
                             MCBSP_RCR_RWDLEN2_16BIT -  16 bits
                             MCBSP_RCR_RWDLEN2_20BIT -  20 bits
                             MCBSP_RCR_RWDLEN2_24BIT -  24 bits
                             MCBSP_RCR_RWDLEN2_32BIT -  32 bits                */

		MCBSP_RCR_RCOMPAND_MSB,/* Receive companding mode (RCOMPAND)
                             MCBSP_RCR_RCOMPAND_MSB  - No companding.Data
                                  transfer starts with MSB first.
                             MCBSP_RCR_RCOMPAND_8BITLSB - No companding,
                                  8-bit data. Transfer starts with LSB first.
                                  Applicable to 8-bit data or 32-bit data in
                                  data reversal mode.
                             MCBSP_RCR_RCOMPAND_ULAW -  Compand using m-law for
                                  receive data. Applicable to 8-bit data only
                             MCBSP_RCR_RCOMPAND_ALAW - Compand using A-law for
                                  receive data. Applicable to 8-bit data only  */

		MCBSP_RCR_RFIG_YES, /* Receive frame ignore(RFIG)
                             MCBSP_RCR_RFIG_NO  - Unexpected receive frame
                                  synchronization pulses restart the transfer.
                             MCBSP_RCR_RFIG_YES - Unexpected receive frame
                                  synchronization pulses are ignored.          */

		MCBSP_RCR_RDATDLY_1BIT,/* Receive data delay(RDATDLY)
                             MCBSP_RCR_RDATDLY_0BIT - 0 bit data delay
                             MCBSP_RCR_RDATDLY_1BIT - 1 bit data delay
                             MCBSP_RCR_RDATDLY_2BIT - 2 bit data delay         */

		MCBSP_RCR_RFRLEN1_OF(0), /* Receive frame length in phase 1(RFRLEN1)
                             000 0000b: 1 word per phase
                             000 0001b: 2 words per phase
                             . . . . . . . . . . . .
                             111 1111b: 128 words per phase                    */

		MCBSP_RCR_RWDLEN1_24BIT,/* Receive element length in phase 1(RWDLEN1)

                             MCBSP_RCR_RWDLEN1_8BIT  -  8 bits
                             MCBSP_RCR_RWDLEN1_12BIT - 12 bits
                             MCBSP_RCR_RWDLEN1_16BIT - 16 bits
                             MCBSP_RCR_RWDLEN1_20BIT - 20 bits
                             MCBSP_RCR_RWDLEN1_24BIT - 24 bits
                             MCBSP_RCR_RWDLEN1_32BIT - 32 bits                 */

		MCBSP_RCR_RWDREVRS_DISABLE/* Receive 32-bit bit reversal feature.(RWDREVRS)
                             MCBSP_RCR_RWDREVRS_DISABLE -32 bit reversal disabled
                             MCBSP_RCR_RWDREVRS_ENABLE  -32 bit reversal enabled.
                                  32-bit data is received LSB first.
                                  RWDLEN should be set for 32-bit operation.
                                  RCOMPAND should be set to 01b else operation
                                  is undefined.                                */
 ),

 MCBSP_XCR_RMK(

    MCBSP_XCR_XPHASE_SINGLE,/* Transmit phases
                             MCBSP_XCR_XPHASE_SINGLE - Single phase frame
                             MCBSP_XCR_XPHASE_DUAL   - Dual phase frame       */

    MCBSP_XCR_XFRLEN2_OF(0),/* Transmit frame length in phase 2(XFRLEN2)
                             000 0000b: 1 word per phase
                             000 0001b: 2 words per phase
                             . . . . . . . . . . . .
                             111 1111b: 128 words per phase                   */

    MCBSP_XCR_XWDLEN2_24BIT, /*  Transmit element length in phase 2
                             MCBSP_XCR_XWDLEN2_8BIT    -  8  bits
                             MCBSP_XCR_XWDLEN2_12BIT   -  12 bits
                             MCBSP_XCR_XWDLEN2_16BIT   -  16 bits
                             MCBSP_XCR_XWDLEN2_20BIT   -  20 bits
                             MCBSP_XCR_XWDLEN2_24BIT   -  24 bits
                             MCBSP_XCR_XWDLEN2_32BIT   -  32 bits             */

    MCBSP_XCR_XCOMPAND_MSB, /* Transmit companding mode(XCOMPAND)
                             MCBSP_XCR_XCOMPAND_MSB      - No companding. Data
                                  transfer starts with MSB first.
                             MCBSP_XCR_XCOMPAND_8BITLSB  - No companding, 8-bit
                                   data. Transfer starts with LSB first.
                                   Applicable to 8-bit data,or 32-bit data in
                                   data reversal mode.
                             MCBSP_XCR_XCOMPAND_ULAW     - Compand using m-law
                                   for receive data.Applicable to 8-bit data only.
                             MCBSP_XCR_XCOMPAND_ALAW  - Compand using A-law for
                                   receive data.Applicable to 8-bit data only. */

    MCBSP_XCR_XFIG_YES, /* Transmit frame ignore(XFIG)
                             MCBSP_XCR_XFIG_NO  - Unexpected transmit frame
                                   synchronization pulses restart the transfer.
                             MCBSP_XCR_XFIG_YES - Unexpected transmit frame
                                   synchronization pulses are ignored.         */

    MCBSP_XCR_XDATDLY_1BIT, /*  Transmit data delay(XDATDLY)
                             MCBSP_XCR_XDATDLY_0BIT   - 0 bit data delay
                             MCBSP_XCR_XDATDLY_1BIT   - 1 bit data delay
                             MCBSP_XCR_XDATDLY_2BIT   - 2 bit data delay       */

    MCBSP_XCR_XFRLEN1_OF(0), /* Transmit frame length in phase 1(XFRLEN1)
                             000 0000b: 1 word per phase
                             000 0001b: 2 words per phase
                             . . . . . . . . . . . .
                             111 1111b: 128 words per phase                    */

    MCBSP_XCR_XWDLEN1_24BIT, /* Transmit element length in phase 1(XWDLEN1)
                             MCBSP_XCR_XWDLEN1_8BIT  -   8 bits
                             MCBSP_XCR_XWDLEN1_12BIT -  12 bits
                             MCBSP_XCR_XWDLEN1_16BIT -  16 bits
                             MCBSP_XCR_XWDLEN1_20BIT -  20 bits
                             MCBSP_XCR_XWDLEN1_24BIT -  24 bits
                             MCBSP_XCR_XWDLEN1_32BIT -  32 bits                */

   MCBSP_XCR_XWDREVRS_DISABLE /* Transmit 32-bit bit reversal feature
                              MCBSP_XCR_XWDREVRS_DISABLE - 32-bit reversal
                                   disabled.
                              MCBSP_XCR_XWDREVRS_ENABLE  - 32-bit reversal
                                   enabled. 32-bit data is transmitted LSB first.
                                   XWDLEN should be set for 32-bit operation.
                                   XCOMPAND should be set to 01b; else operation
                                   is undefined.                               */

  ), /* Transmitter Control Reg. (XCR) */

MCBSP_SRGR_RMK( /* Sample�Rate Generator Reg. (SRGR) */

		MCBSP_SRGR_GSYNC_FREE,/* Sample rate generator clock synchronization(GSYNC).
                           MCBSP_SRGR_GSYNC_FREE - The sample rate generator
                                clock CLKG) is free running.
                           MCBSP_SRGR_GSYNC_SYNC - (CLKG) is running but is
                                 resynchronized, and the frame sync signal
                                 (FSG)is generated only after the receive
                                 frame synchronization signal(FSR)is detected.
                                 Also,the frame period (FPER) is a don�t care
                                 because the period is dictated by the external
                                 frame sync pulse.                            */

		MCBSP_SRGR_CLKSP_RISING,/* CLKS polarity clock edge select(CLKSP)
                            MCBSP_SRGR_CLKSP_RISING  - The rising edge of CLKS
                                 generates CLKG and FSG.
                            MCBSP_SRGR_CLKSP_FALLING - The falling edge of CLKS
                                 generates CLKG and FSG.                      */
		MCBSP_SRGR_CLKSM_INTERNAL,/* MCBSP sample rate generator clock mode(CLKSM)
                            MCBSP_SRGR_CLKSM_CLKS   - The sample rate generator
                                 clock is derived from CLKS.
                            MCBSP_SRGR_CLKSM_INTERNAL - (Default value) The
                                 sample rate generator clock is derived from
                                 the internal clock source.                   */

		MCBSP_SRGR_FSGM_FSG,/*Sample rate generator transmit frame synchronization
                             mode.(FSGM)
                            MCBSP_SRGR_FSGM_DXR2XSR  - The transmit frame sync
                                 signal (FSX) is generated on every DXR to XSR
                                 copy.
                            MCBSP_SRGR_FSGM_FSG      - The transmit frame sync
                                 signal is driven by the sample rate generator
                                 frame sync signal, FSG.                      */

		MCBSP_SRGR_FPER_OF(24),/* Frame period(FPER)
                            Valid values: 0 to 4095                           */

		MCBSP_SRGR_FWID_OF(0),/* Frame width(FWID)
                            Valid values: 0 to 255                            */

		MCBSP_SRGR_CLKGDV_OF(5)/* Sample rate generator clock divider(CLKGDV)
                            Valid values: 0 to 255                            */

),
0x00000000, /* Multichannel Control Reg. (MCR) */
0x00000000, /* Receiver Channel Enable(RCER) */
0x00000000, /* Transmitter Channel Enable(XCER) */
MCBSP_PCR_RMK(

   MCBSP_PCR_XIOEN_SP, /* Transmitter in general-purpose I/O mode - only when
                          XRST = 0 in SPCR - (XIOEN)
                          MCBSP_PCR_XIOEN_SP    -  CLKS pin is not a general
                               purpose input. DX pin is not a general purpose
                               output.FSX and CLKX are not general-purpose I/Os.
                          MCBSP_PCR_XIOEN_GPIO  -  CLKS pin is a general-purpose
                               input. DX pin is a general-purpose output.
                               FSX and CLKX are general-purpose I/Os. These
                               serial port pins do not perform serial port
                               operation.                                     */
   MCBSP_PCR_RIOEN_SP, /* Receiver in general-purpose I/O mode - only when
                          RRST = 0 in SPCR -(RIOEN)
                          MCBSP_PCR_RIOEN_SP    - DR and CLKS pins are not
                               general-purpose inputs. FSR and CLKR are not
                               general-purpose I/Os and perform serial port
                               operation.
                          MCBSP_PCR_RIOEN_GPIO  - DR and CLKS pins are
                               general-purpose inputs. FSR and CLKR are
                               general-purpose I/Os. These serial port pins do
                               not perform serial port operation.            */
   MCBSP_PCR_FSXM_INTERNAL, /* Transmit frame synchronization mode(FSXM)
                            MCBSP_PCR_FSXM_EXTERNAL - Frame synchronization
                                 signal is provided by an external source. FSX
                                 is an input pin.
                            MCBSP_PCR_FSXM_INTERNAL - Frame synchronization
                                 generation is determined by the sample rate
                                 generator frame synchronization mode bit FSGM
                                 in the SRGR.                                */

   MCBSP_PCR_FSRM_INTERNAL, /* Receive frame synchronization mode (FSRM)
                            MCBSP_PCR_FSRM_EXTERNAL  - Frame synchronization
                                 signals are generated by an external device.
                                 FSR is an input pin.
                            MCBSP_PCR_FSRM_INTERNAL  - Frame synchronization
                                 signals are generated internally by the sample
                                 rate generator. FSR is an output pin except
                                 when GSYNC = 1 in SRGR.                     */

   MCBSP_PCR_CLKXM_OUTPUT, /* Transmitter clock mode (CLKXM)
                            MCBSP_PCR_CLKXM_INPUT    -  Transmitter clock is
                                 driven by an external clock with CLKX as an
                                 input pin.
                            MCBSP_PCR_CLKXM_OUTPUT   - CLKX is an output pin
                                  and is driven by the internal sample rate
                                  generator.

                            During SPI mode :
                            MCBSP_PCR_CLKXM_INPUT    -  McBSP is a slave and
                                 (CLKX) is driven by the SPI master in the
                                  system. CLKR is internally driven by CLKX.
                            MCBSP_PCR_CLKXM_OUTPUT   - McBSP is a master and
                                 generates the transmitter clock (CLKX) to
                                 drive its receiver clock (CLKR) and the shift
                                 clock of the SPI-compliant slaves in the
                                 system.                                     */
   MCBSP_PCR_CLKRM_INPUT, /* Receiver clock mode (CLKRM)

                             Case 1: Digital loopback mode not set in SPCR

                             MCBSP_PCR_CLKRM_INPUT - Receive clock (CLKR) is
                                  an input driven by an external clock.

                             MCBSP_PCR_CLKRM_OUTPUT -  CLKR is an output pin
                                  and is driven by the sample rate generator.

                             Case 2: Digital loopback mode set  in SPCR

                             MCBSP_PCR_CLKRM_INPUT - Receive clock  is driven
                                  by the transmit clock (CLKX), which is based
                                  on the CLKXM bit in PCR. CLKR is in high
                                  impedance.
                             MCBSP_PCR_CLKRM_INPUT - CLKR is an output pin and
                                  is driven by the transmit clock. The transmit
                                  clock is derived from CLKXM bit in the PCR.*/


   MCBSP_PCR_CLKSSTAT_0, /*  CLKS pin status(CLKSSTAT)
                             MCBSP_PCR_CLKSSTAT_0
                             MCBSP_PCR_CLKSSTAT_1                            */

   MCBSP_PCR_DXSTAT_0,   /*  DX pin status(DXSTAT)
                             MCBSP_PCR_DXSTAT_0
                             MCBSP_PCR_DXSTAT_1                              */

   MCBSP_PCR_FSXP_ACTIVEHIGH, /* Transmit frame synchronization polarity(FSXP)
                             MCBSP_PCR_FSXP_ACTIVEHIGH - Frame synchronization
                                      pulse FSX is active high
                             MCBSP_PCR_FSXP_ACTIVELOW  - Frame synchronization
                                  pulse FSX is active low                    */
   MCBSP_PCR_FSRP_ACTIVEHIGH, /* Receive frame synchronization polarity(FSRP)
                             MCBSP_PCR_FSRP_ACTIVEHIGH - Frame synchronization
                                  pulse FSR is active high
                             MCBSP_PCR_FSRP_ACTIVELOW  - Frame synchronization
                                  pulse FSR is active low                    */
   MCBSP_PCR_CLKXP_RISING, /* Transmit clock polarity(CLKXP)
                             MCBSP_PCR_CLKXP_RISING - Transmit data driven on
                                  rising edge of CLKX
                             MCBSP_PCR_CLKXP_FALLING - Transmit data driven on
                                   falling edge of CLKX                      */
   MCBSP_PCR_CLKRP_RISING /* Receive clock polarity(CLKRP)
                             MCBSP_PCR_CLKRP_FALLING - Receive data sampled on
                                  falling edge of CLKR
                             MCBSP_PCR_CLKRP_RISING - Receive data sampled on
                                   rising edge of CLKR                       */
 ) /* Pin Control Reg. (PCR) */
};


#endif
