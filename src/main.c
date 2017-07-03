#include <stdio.h>
#include <stdlib.h>
#include <magnetometer.h>
#include "libcsv/csv.h"
#include "AccelAndMagnet/magnetometer.h"

int gi = 0;
uvector_t points[100];
void cb1(void *s, size_t i, void **data)
{
    float temp = atof(s);
    ((uvector_t *)*data)->arr[gi] = F16(temp);
    gi++;
}

void cb2(int c, void **data)
{
    (*(uvector_t **)data)++;
    gi = 0;
}

int main(int argc, char* argv[])
{
    struct csv_parser p;
    uvector_t *pPoints = &points[0];
    unsigned char options = 0;
    csv_init(&p, options);
    csv_set_delim(&p, ',');
    char buffer[1024];
    size_t bytes_read;
    FILE *fstream = fopen("data.txt","r");
    while ((bytes_read=fread(buffer, 1, 1024, fstream)) > 0) {
        if (csv_parse(&p, buffer, bytes_read, &cb1, &cb2, (void **)&pPoints) != bytes_read) {
            fprintf(stderr, "Error while parsing file: %s\n", csv_strerror(csv_error(&p)));
        }
    }

    csv_fini(&p, &cb1, &cb2, (void **) &pPoints);
    return 0;

}