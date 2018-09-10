#ifndef __PSPMTSTRUCT_HPP__
#define __PSPMTSTRUCT_HPP__

struct PspmtEvent {
        double x_position = -1;
        double y_position = -1;
        int x_pixel = -1;
        int y_pixel = -1;
        int pixel_num = -1;
        double event_time = -1;
        bool implant = false;
        bool decay = false;
        double low_dynode = 0;
        double hi_dynode = 0;
        double low_dynode_mult = 0;
        double hi_dynode_mult = 0;

    };
#endif
